#include "namespace.h"
#include "dataIO.h"
#include "gibbsSampler.h"
#include "utilities.h"

using namespace std;

GibbsSampler :: GibbsSampler(DataIO &dataIOObj)
{
	baseAlpha = 0.01;
	baseBeta = 1;
	logBase = 1.5;
	logLikelihood = 0;

	initializeForSampler(dataIOObj);

	updateNodeNodeCountMap(dataIOObj);
	
	initializeUserUserInfluence(dataIOObj);
	initializeAvgProbabilityVectors(dataIOObj);
	initializeAvgTopicProbabilityVectors(dataIOObj);

	// So that its not just the prior probabilites...
	// So initializing/or sampling once to have an effect of the initial counts...
	sampleInfluenceAssignment(dataIOObj);
	// update the user base rates...
	updateUserBaseRates(dataIOObj);
}


int GibbsSampler :: IterativeSampler(DataIO &dataIOObj, ui ITERATIONS, int BURN_IN)
{
	for(ui ITE = 0; ITE < ITERATIONS; ITE++)
	{
		logLikelihood = 0;

		std::chrono::high_resolution_clock::time_point t1, t2;
		auto duration = 0;

		// we are not dealing with topic assignment... In this case the topic assignments are known...
		 
		t1 = std::chrono::high_resolution_clock::now();
		sampleTopicAssignment(dataIOObj, ITE, BURN_IN);
		t2 = std::chrono::high_resolution_clock::now();
	
		duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
		cout << "sampled topics for iteration -- "<< ITE << "--  Time taken -- " << duration << endl;


		t1 = std::chrono::high_resolution_clock::now();
		sampleParentAssignment(dataIOObj, ITE, BURN_IN);
		t2 = std::chrono::high_resolution_clock::now();
	
		duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
		cout << "sampled parents for iteration -- "<< ITE << "--  Time taken -- " << duration << endl;


		// updating the childEventsMap... as it might change after each iteration once we incorporate parent sampling
		t1 = std::chrono::high_resolution_clock::now();

		dataIOObj.childEventsMap.clear();
		for(unsigned int i = 0; i < dataIOObj.newSyntheticEvents.size(); i++)
		{
			li eventCurrParent = dataIOObj.newSyntheticEvents[i][2];
			if(eventCurrParent != -1)
			{
				dataIOObj.childEventsMap[eventCurrParent].push_back(i);
			}
		}

		updateNodeNodeCountMap(dataIOObj);

		t2 = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
		cout << "Updating Child Events Map and Node-Node Count after -- "<< ITE << " Iteration --  Time  taken -- " << duration << endl;



		cout << setprecision(10) << "Log Likelihood = " << logLikelihood << "\n";

		// llFile << setprecision(10) << logLikelihood << "\n";

		cout << "ITE -- " << ITE << endl;
		
		if(ITE % 100 == 0)
		{
			cout << "ITE -- " << ITE << endl;
		}

	}

	return 0;
}


// vector< vector <li> > GibbsSampler :: initializeForSampler(vector< vector <li> > allEvents, int flag)
int GibbsSampler :: initializeForSampler(DataIO &dataIOObj)
{
	// Random topic initialization...
	// Parent Initialization is time based parent...
	vector< vector <li> > localNewSyntheticEvents;
	// ll uNode, vNode;
	// int eventIndex = 0;
	
	for(ui i = 0; i < dataIOObj.allEvents.size(); i++)
	{
		vector<li> tempEvent;
		tempEvent.push_back(i);
		tempEvent.push_back(dataIOObj.allEvents[i][1]);
		// tempEvent.push_back(-2);											// some arbit value.... saying that its hidden

		int currentEventNode = dataIOObj.allEvents[i][1];

		if(dataIOObj.allPossibleParentEvents[i].size() > 0)
		{
			// initialize to nearest in time parent...
			int possParentEventId = dataIOObj.allPossibleParentEvents[i][0];

			tempEvent.push_back(possParentEventId);				// some initialization at the start.... i.e. nearest in time from neighbor...

			// populate the childevents for the parent events...
			dataIOObj.childEventsMap[possParentEventId].push_back(i);

			int possParentEventNode = localNewSyntheticEvents[possParentEventId][1];
			// checking or updating the map that keeps track of the updated edges...
			dataIOObj.nodeNodeCombUpdateInfluence[possParentEventNode].push_back(currentEventNode);
			// getting the Nuv value... will be updated after every iteration...
			dataIOObj.nodeNodeCount[possParentEventNode][currentEventNode]++;
		}
		else
		{
			// if no nearest in time parent, assign no parent...
			tempEvent.push_back(-1);
		}
		// assigning random topic...
		tempEvent.push_back(i % dataIOObj.numTopics);										// for some arbit initialization...

		localNewSyntheticEvents.push_back(tempEvent);
		
		// populate the node and eventIndex map... will be useful while searching for possible parent events...
		dataIOObj.nodeEventsMap[(ui)currentEventNode].push_back(i);
		// the events here must be already in sorted order, because of the way they are added....

		// i = eventIndex...
		// if (levelInfo[i] != maxLevel)
		
		
		if (i <= 150000)
		{
			dataIOObj.nodeEventsCountMap[(ui)currentEventNode]++;
		}
	}

	dataIOObj.newSyntheticEvents = localNewSyntheticEvents;

	updateCountMatrices(dataIOObj, localNewSyntheticEvents);

	return 0;
}

int GibbsSampler :: updateCountMatrices(DataIO &dataIOObj, vector< vector <li> > &localNewSyntheticEvents)
{
	// initialization of the NTT Matrix... because of some initialization of parents...
	int sumAll = 0;
	for(ui i = 0; i < localNewSyntheticEvents.size(); i++)
	{
		int eParent = localNewSyntheticEvents[i][2];
		int eTopic = localNewSyntheticEvents[i][3];

		if(eParent > -1)
		{
			int eParentTopic = localNewSyntheticEvents[eParent][3];
			dataIOObj.NTTCountVec[eParentTopic][eTopic] += 1;
			dataIOObj.NTTSumTopicsVec[eParentTopic] += 1;
			sumAll++;
		}
		else
		{
			int eNode = localNewSyntheticEvents[i][1];
			dataIOObj.NUTCountVec[eNode][eTopic] += 1;
			dataIOObj.NUTSumTopicsVec[eNode] += 1;
		}
	}

	cout << "Number of Edges = " << sumAll << "\n";

	// populate the topic-word matrix as well as per the topic assignment...
	for(unsigned int i = 0; i < dataIOObj.allDocs.size(); i++)
	{
		// assignedSampleTopics.push_back(-1);
		int initTopic = localNewSyntheticEvents[i][3];;

		for(ui j = 0; j < dataIOObj.allDocs[i].size(); j++)
		{
			int word = dataIOObj.allDocs[i][j];
			dataIOObj.NTWCountVec[initTopic][word]++;
		}

		dataIOObj.NTWSumWordsVec[initTopic] += dataIOObj.allDocs[i].size();
	}

	return 0;
}


int GibbsSampler :: updateNodeNodeCountMap(DataIO &dataIOObj)
{
	li eventNode, parentEvent, parentNode;
	ui Nuv, Nu;
	double alphaDash, betaDash;
	ui uNode, vNode;

	dataIOObj.nodeNodeCount.clear();
	dataIOObj.actualEdgesNuSum.clear();
	// triggeredEdges.clear();
	dataIOObj.groupTransactionsSum.clear();

	for(unsigned int i = 0; i < dataIOObj.newSyntheticEvents.size(); i++)
	{
		eventNode = dataIOObj.newSyntheticEvents[i][1];
		parentEvent = dataIOObj.newSyntheticEvents[i][2];
		if(parentEvent >= 0)
		{
			parentNode = dataIOObj.newSyntheticEvents[parentEvent][1];
			dataIOObj.nodeNodeCount[parentNode][eventNode]++;
		}
	}

	// unordered_map<ui, unordered_map<ui, ui> >::iterator nodeNodeCountIterator;
	map<ui, map<ui, ui> >::iterator nodeNodeCountIterator;

	// this is an overhead.. as we are calculating the values in the values while estimating the influence...
	for(nodeNodeCountIterator = dataIOObj.nodeNodeCount.begin(); nodeNodeCountIterator != dataIOObj.nodeNodeCount.end(); nodeNodeCountIterator++)
	{
		uNode = nodeNodeCountIterator->first;
		Nu = dataIOObj.nodeEventsCountMap[uNode];

		// unordered_map<ui, ui> nodeCountMap = nodeNodeCountIterator->second;
		map<ui, ui> nodeCountMap = nodeNodeCountIterator->second;
		// unordered_map<ui, ui>::iterator nodeCountMapIterator;
		map<ui, ui>::iterator nodeCountMapIterator;

		for(nodeCountMapIterator = nodeCountMap.begin(); nodeCountMapIterator != nodeCountMap.end(); nodeCountMapIterator++)
		{
			vNode = nodeCountMapIterator->first;
			Nuv = nodeCountMapIterator->second;

			int outDeg = dataIOObj.outDegreeMap[uNode];
			int inDeg = dataIOObj.inDegreeMap[vNode];

			int outDegGID = floor(log2(outDeg)/log2(logBase));
			int inDegGID = floor(log2(inDeg)/log2(logBase));

			string gid = to_string(outDegGID).append("_").append(to_string(inDegGID));

			dataIOObj.groupTransactionsSum[gid] += Nuv;

			dataIOObj.actualEdgesNuSum[gid] += Nu;
            
            // sumPoissonLambda[gid] += (Nuv*1.0)/Nu;
            // countSumPoissonLambda[gid] += 1;

			// triggeredEdges[gid] += 1;
		}
	}

	return 0;
}


int GibbsSampler :: initializeUserUserInfluence(DataIO &dataIOObj)
{
	cout << "Initializing user user influence to a prior...\n";

	ui count = 0;
	ui uNode, vNode;
	// int numTimesAddedNu = 0;

	for(ui i = 0; i < dataIOObj.followersMap.size(); i++)
	{
		vector <ui> followers = dataIOObj.followersMap[i];

		uNode = i;
		
		unordered_map <ui, double> tempNodeInf;

		int sourceOutDeg = followers.size();
		int sourceGid;

		if(sourceOutDeg > 0)
		{
			sourceGid = (int) floor(log2(sourceOutDeg)/log2(logBase));
			for(ui j = 0; j < followers.size(); j++)
			{
				vNode = followers[j];
				double scaleParam = 1/(baseBeta);

				// tempNodeInf[vNode] = getSampleFromGamma(baseAlpha, scaleParam);
				tempNodeInf[vNode] = baseAlpha * scaleParam;

				int destInDeg = dataIOObj.inDegreeMap[vNode];
				int destGid = (int) floor(log2(destInDeg)/log2(logBase));

				string gid = to_string(sourceGid).append("_").append(to_string(destGid));

				// one time calcualtion...
				// as the number of events for each node do not change over iterations....
				// also, the contribution is taken only from the first 15% of the events... 
				dataIOObj.groupSourceTransactionsSum[gid] += dataIOObj.nodeEventsCountMap[uNode];
				
				// edgesInGroup[gid] += 1;
				// numTimesAddedNuToGroup[gid] += 1;
				// numTimesAddedNu += 1;

				// cout << nodeEventsCountMap[uNode] << endl;
			}
			dataIOObj.userUserInfluence[uNode] = tempNodeInf;
		}

		tempNodeInf.clear();

		count++;

		if(count % 100000 == 0)
			cout << "Initialized for " <<  count << " users...\n";
	}
    dataIOObj.followersMap.clear();
	cout << "Initialized user user influence to a prior...\n";
	// cout << "Added Nu number of times = " << numTimesAddedNu << endl;
	return 0;
}


int GibbsSampler ::  initializeAvgProbabilityVectors(DataIO &dataIOObj)
{
	for(ui i = 0; i < dataIOObj.allPossibleParentEvents.size(); i++)
	{
		int probVecSize = dataIOObj.allPossibleParentEvents[i].size();
		vector<double> initVec(probVecSize + 1, 0.0);
		dataIOObj.avgProbParForAllEvents.push_back(initVec);
		initVec.clear();
	}

	return 0;
}

int GibbsSampler :: initializeAvgTopicProbabilityVectors(DataIO &dataIOObj)
{
	for(ui i = 0; i < dataIOObj.allDocs.size(); i++)
	{
		vector <double> initVec(dataIOObj.numTopics, 0.0);
		dataIOObj.avgTopicProbVector.push_back(initVec);
		initVec.clear();
	}
	return 0;
}





///////////////////////////////
/////// User Base Rates  //////
///////////////////////////////
int GibbsSampler :: updateUserBaseRates(DataIO &dataIOObj)
{
	int sponCount = 0;
	int minTime = dataIOObj.eventIndexTimestamps[0];
	int maxTime = dataIOObj.eventIndexTimestamps[dataIOObj.eventIndexTimestamps.size() - 1];

	int totalDataObservedTime = maxTime - minTime;

	map<int, int> eachNodeSponCount;


	map<int, bool> distinctNodes;

	for(ui i = 0; i < dataIOObj.newSyntheticEvents.size(); i++)
	{
		if(dataIOObj.newSyntheticEvents[i][2] == -1)
		{
			sponCount++;
			distinctNodes[dataIOObj.newSyntheticEvents[i][1]] = 1;
			eachNodeSponCount[dataIOObj.newSyntheticEvents[i][1]]++;
		}
	}

	// defaultMuVal = (sponCount * 1.0) / ((maxTime - minTime) * distinctNodes.size());

	double averageMuVal = 0;

	double minMuVal = 1000000;
	double maxMuVal = 0;
	map<int, int>::iterator eachNodeSponCountIt;	

	for(eachNodeSponCountIt = eachNodeSponCount.begin(); eachNodeSponCountIt != eachNodeSponCount.end(); eachNodeSponCountIt++)
	{
		int nodeId = eachNodeSponCountIt->first;
		int nodeSponCount = eachNodeSponCountIt->second;

		dataIOObj.userBaseRateMap[nodeId] = (nodeSponCount * 1.0) / totalDataObservedTime;

		averageMuVal +=  dataIOObj.userBaseRateMap[nodeId];

		if(dataIOObj.userBaseRateMap[nodeId] > maxMuVal)
		{
			maxMuVal = dataIOObj.userBaseRateMap[nodeId];
		}

		if(dataIOObj.userBaseRateMap[nodeId] < minMuVal)
		{
			minMuVal = dataIOObj.userBaseRateMap[nodeId];
		}
	}

	dataIOObj.defaultMuVal = averageMuVal / eachNodeSponCount.size();
	
	cout << "sponCount = " << sponCount << " minMuVal = " << minMuVal << " maxMuVal = " << maxMuVal << " MLE MuVal = " << dataIOObj.defaultMuVal << endl;

	return 0;
}


///////////////////////////////
/////// SAMPLE INFLUENCE //////
///////////////////////////////
/*
int GibbsSampler :: sampleInfluenceAssignment(DataIO &dataIOObj)
{
	cout << "Sample User User Influence\n";
	ui Nuv, Nu;
	double alphaDash, betaDash;
	ui uNode, vNode;

	unordered_map<ui, vector<ui> >::iterator nodeNodeUpdateIterator;

	double avgVal = 0;
	int count = 0;

	for(nodeNodeUpdateIterator = dataIOObj.nodeNodeCombUpdateInfluence.begin(); nodeNodeUpdateIterator != dataIOObj.nodeNodeCombUpdateInfluence.end(); nodeNodeUpdateIterator++)
	{
		uNode = nodeNodeUpdateIterator->first;
		vector<ui> updateForNodes = nodeNodeUpdateIterator->second;

		Nu = nodeEventsCountMap[uNode];

		for(ui i = 0; i < updateForNodes.size(); i++)
		{
			vNode = updateForNodes[i];
			Nuv = dataIOObj.nodeNodeCount[uNode][vNode];

			alphaDash = Nuv + dataIOObj.baseAlpha;
			// betaDash = 1/(Nu + dataIOObj.baseBeta);
			// double scaleParam = 1/(Nu + rateParam);
			double scaleParam = 1/(Nu + dataIOObj.baseBeta);
			// betaDash = 1.0/Nu;                 			//remove

			// userUserInfluence[uNode][vNode] = getSampleFromGamma(alphaDash, scaleParam);
			userUserInfluence[uNode][vNode] = alphaDash * scaleParam;

			avgVal += dataIOObj.userUserInfluence[uNode][vNode];
			count++;
		}
	}

	cout << "Avg Wuv = " << avgVal / count << endl;

	return 0;
}
*/

int GibbsSampler :: sampleInfluenceAssignment(DataIO &dataIOObj)
{
	cout << "Sample User User Influence\n";
	ui Nuv, Nu;
	double alphaDash, betaDash;
	ui uNode, vNode;

	unordered_map<ui, vector<ui> >::iterator nodeNodeUpdateIterator;

	double avgVal = 0;
	int count = 0;

	// unordered_map<ui, unordered_map<ui, ui> >::iterator nodeNodeCountIterator;
	map<ui, map<ui, ui> >::iterator nodeNodeCountIterator;

	for(nodeNodeCountIterator = dataIOObj.nodeNodeCount.begin(); nodeNodeCountIterator != dataIOObj.nodeNodeCount.end(); nodeNodeCountIterator++)
	{
		uNode = nodeNodeCountIterator->first;
		int outDeg = dataIOObj.outDegreeMap[uNode];

		// vector<ui> updateForNodes = nodeNodeCountIterator->second;
		// unordered_map<ui, ui> nodeCountMap = nodeNodeCountIterator->second;
		map<ui, ui> nodeCountMap = nodeNodeCountIterator->second;		
		// unordered_map<ui, ui>::iterator nodeCountMapIterator;		
		map<ui, ui>::iterator nodeCountMapIterator;		

		for(nodeCountMapIterator = nodeCountMap.begin(); nodeCountMapIterator != nodeCountMap.end(); nodeCountMapIterator++)
		{
			vNode = nodeCountMapIterator->first;
			// Nuv = nodeCountMapIterator->second;

			int inDeg = dataIOObj.inDegreeMap[vNode];

			int outDegGID = floor(log2(outDeg)/log2(logBase));
			int inDegGID = floor(log2(inDeg)/log2(logBase));

			string gid = to_string(outDegGID).append("_").append(to_string(inDegGID));

			// default uinf 
			double uinf = 0.01;

			// if the denominator is not zero... 
			// This can happen because we are taking the contribution to the denominator 
			// of Wuv only from the first 15% of the events... 
			// if an event which is not a part of first 15% of the events
			// then the Nu count for this event would be zero... 
			// and might have zero for the group... 
			// the chances are lesser here because of group formations...
			if(dataIOObj.groupSourceTransactionsSum[gid] != 0)
			{
				uinf = (dataIOObj.groupTransactionsSum[gid] + baseAlpha) * 1.0 / (dataIOObj.groupSourceTransactionsSum[gid] + baseBeta);
			}

			dataIOObj.groupWuv[gid] = uinf;

			dataIOObj.userUserInfluence[uNode][vNode] = uinf;

			avgVal += uinf;
			count += 1;
		}
	}

	cout << "Avg Wuv = " << avgVal / count << endl;

	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////
///////  SAMPLE TOPIC    //////
///////////////////////////////

int GibbsSampler :: sampleTopicAssignment(DataIO &dataIOObj, int ITE, int BURN_IN)
{
	cout << "Sample Topic Assignment\n";

	for(unsigned int i = 0; i < dataIOObj.newSyntheticEvents.size(); i++)	
	{
		li eventIndex = dataIOObj.newSyntheticEvents[i][0];
		// li eventIndex = i;
		li eventNode = dataIOObj.newSyntheticEvents[i][1];
		li eventParent = dataIOObj.newSyntheticEvents[i][2];			// should be -2 for the first iteration...
		li eventTopic = dataIOObj.newSyntheticEvents[i][3];				// should be -2 for the first iteration...	

		vector<ui> doc = dataIOObj.allDocs[i];

		li assignedTopic;
		
		decreamentCountFromMatrices(dataIOObj, eventIndex, eventNode, eventParent, eventTopic, doc, 1);
		assignedTopic =  getSampledTopicAssignment(dataIOObj, eventIndex, eventNode, eventParent, doc, ITE, BURN_IN);
	
		increamentCountToMatrices(dataIOObj, eventIndex, eventNode, eventParent, assignedTopic, doc, 1);

		dataIOObj.newSyntheticEvents[eventIndex][3] = assignedTopic;

		// getting every 10th iteration topic assigned, i.e. the topic assigned at each Markov Chain... 
		if(ITE > BURN_IN && ITE%10 == 0)
		{
			dataIOObj.topicDistEvery10thIter[eventIndex].push_back(assignedTopic);
		}
	}

	printTopicTopicCount(dataIOObj);

	return 0;
}


li GibbsSampler :: getSampledTopicAssignment(DataIO &dataIOObj, li eventIndex, li eventNode, li eventParent, vector<ui> doc, int ITE, int BURN_IN)
{
	li assignedTopic = -1;

	vector <double> calculatedProbVec (dataIOObj.numTopics, 0.0);
	// we need the hist of topics assigned to the child events... For the second/middle term specifically...
	unordered_map <int, ui> childEventTopicsHist = getHistOfTopicsOverChildEvents(dataIOObj, eventIndex);

	// calculating probability for each possible topic
	for(li topic = 0; topic < dataIOObj.numTopics; topic++)
	{
		// cout << "Topic iterator at - " << topic << endl;
		double firstTerm = 0.0;
		double middleTerm = 0.0;
		double thirdTerm = 0.0;

		// cout << "eventIndex -- " << eventIndex << endl;
		// first term -- related to the parents topic... or else to the node Topic preference in case of no parent...
		firstTerm = getFirstTermOfTopicAssignmentCondProb(dataIOObj, eventNode, eventParent, topic);
		
		if(childEventTopicsHist.size() > 0)
		{
			middleTerm = getMiddleTermOfTopicAssignmentCondProb(dataIOObj, topic, childEventTopicsHist, eventIndex, eventParent);
		}
		else
		{
			// if there is no child event... then there is no event to iterate on... 
			// and the probability should be calculated just on the basis of the parent event and words in the event...
			middleTerm = 0;
		}
		// middleTerm = 0;

		// hist of words for each doc do not change over iteration...  
		vector<ui> wordHistVec = dataIOObj.wordHistAllDocsVector[eventIndex];
		thirdTerm = getThirdTermOfTopicAssignmentCondProb(dataIOObj, topic, doc, wordHistVec);
		// thirdTerm = 0;

		calculatedProbVec[topic] = firstTerm + middleTerm + log(thirdTerm);
		// cout << "Topic iterator at done -- " << topic << endl;
	}

	assignedTopic = getSampleFromMultinomial(calculatedProbVec);

	if(ITE > BURN_IN && ITE%10 == 0)
	{
		for(ui i = 0; i < dataIOObj.avgTopicProbVector[eventIndex].size(); i++)
		{
			dataIOObj.avgTopicProbVector[eventIndex][i] += calculatedProbVec[i];
		}
	}

	calculatedProbVec.clear();
	// cout << assignedTopic << "\n";

	return assignedTopic;
}

double GibbsSampler :: getFirstTermOfTopicAssignmentCondProb(DataIO &dataIOObj, li eventNode, li eventParent, li topic)
{
	double finalFirstTerm = 0;

	if(eventParent >= 0)
	{
		li eventParentTopic = dataIOObj.newSyntheticEvents[eventParent][3];
		// finalFirstTerm = hyperBeta[topic] + NTTCountVec[eventParentTopic][topic];
		finalFirstTerm = dataIOObj.hyperBeta[topic] + dataIOObj.NTTCountVec[eventParentTopic][topic];
	}
	else if(eventParent == -1)
	{
		// no parent
		finalFirstTerm = dataIOObj.hyperGamma[topic] + dataIOObj.NUTCountVec[eventNode][topic];
		// cout << eventNode << " " << topic << endl;
	}
	else
	{
		// This would not get executed if there is some parent initialization...
		// if the parent is invalid or not assigned yet...
		// finalFirstTerm = log(hyperGamma[topic]);
		finalFirstTerm = dataIOObj.hyperGamma[topic];
	}

	double logVal = log(finalFirstTerm);
	// cout << logVal << endl;
	return logVal;
}


double GibbsSampler :: getMiddleTermOfTopicAssignmentCondProb(DataIO &dataIOObj, li topic, unordered_map <int, ui> childEventTopicsHist, li eventIndex, li eventParent)
{
	double finalMiddleValue = 0.0;
	double finalMiddleValueNume = 0.0;
	double finalMiddleValueDenom = 0.0;

	li eventParentTopic;

	if(eventParent >= 0)
		eventParentTopic = dataIOObj.newSyntheticEvents[eventParent][3];
	else
		eventParentTopic = -2;				// some invalid value...

	unordered_map <int, ui>::iterator childEventTopicsHistIt;

	// int childEventTopicsHistCount = 0;
	for(childEventTopicsHistIt = childEventTopicsHist.begin(); childEventTopicsHistIt != childEventTopicsHist.end(); childEventTopicsHistIt++)
	{
		// count of the topics in the child events...
		li lPrime = childEventTopicsHistIt->first;

		// if there is some topic assigned to the child events...  
		if(lPrime > -1)
		{
			// For a particular possible topic the base term remains same 
			// over the count of the topics in the neighborhood/child events...
			li lPrimeCountInChilds = childEventTopicsHistIt->second;

			// double baseTerm = betaLPrime + ttCountWithoutChildEvents;
			double baseTerm = dataIOObj.hyperBeta[lPrime] + dataIOObj.NTTCountVec[topic][lPrime];

			double valueForLPrime = 0.0;

			// evaluating sum( log ( baseterm + i) )

			for(li i = 0; i < lPrimeCountInChilds; i++)
			{
				valueForLPrime += log(baseTerm + i);
			}

			finalMiddleValueNume += valueForLPrime;
		}
	}

	// Denominator is to be calculated over the number of child events...
	// the baseterm remains same for all the terms...

	double baseTermDenom = dataIOObj.sumBeta + dataIOObj.NTTSumTopicsVec[topic];

	ui numChildEvents = 0;

	vector <ui> cEvents;

	try
	{
		// get the number of child events from the updated child events map...
		numChildEvents = dataIOObj.childEventsMap[eventIndex].size();
		// remove the count of the child events that have an unassigned topic...
		// this should not matter if there is some initialization... 
		numChildEvents = numChildEvents - childEventTopicsHist[-1];

	}
	catch(exception &e)
	{
		numChildEvents = 0;
	}

	if(numChildEvents > 0)
	{
		for(ui i = 0; i < numChildEvents; i++)
		{
			finalMiddleValueDenom += log(baseTermDenom + i);
		}
	}
	else
	{
		finalMiddleValueDenom = 0;
	}
	// it cannot happen that denominator is zero and numerator is not... 
	if(finalMiddleValueDenom == 0.0 && finalMiddleValueNume != 0.0)
	{
		cout << "Some issue with the code...\n";
		exit(0);
	}

	finalMiddleValue = finalMiddleValueNume - finalMiddleValueDenom;

	return finalMiddleValue;
}

// this is required for calculating the middleterm of the propobability for topic assignment...
unordered_map <int, ui>  GibbsSampler :: getHistOfTopicsOverChildEvents(DataIO &dataIOObj, li eventIndex)
{
	unordered_map <int, ui> topicCounts;

	vector <ui> childEventsList;
	try
	{
		childEventsList =  dataIOObj.childEventsMap.at(eventIndex);
	}
	catch(exception &e)
	{
		return topicCounts;							// returning empty dict/map
	}

	if(childEventsList.size() > 0)
	{
		for(unsigned int i = 0; i < childEventsList.size(); i++)
		{
			ui childEventIndex = childEventsList[i];
			li childEventTopic = dataIOObj.newSyntheticEvents[childEventIndex][3];

			topicCounts[childEventTopic] += 1;
		}
	}

	return topicCounts;
}  


double  GibbsSampler :: getThirdTermOfTopicAssignmentCondProb(DataIO &dataIOObj, li topic, vector<ui> doc, vector<ui> wordHistVec)
{
	// as the documents are short we are not doing log computation here... and calculate log after returning the computed value...
	long double finalThirdValue = 1.0;
	long double finalThirdValueNume = 1.0;
	long double finalThirdValueDenom = 1.0;

	unsigned i = 0;
	while(i < wordHistVec.size())
	{
		ui word = wordHistVec[i];
		i++;
		ui wordCount = wordHistVec[i];

		// as for all the words this is same... we do not really need hyperAlpha right now...
		// double baseterm = alphaWord + NTWCountVec[topic][word];
		double baseterm = dataIOObj.hyperAlpha[word] + dataIOObj.NTWCountVec[topic][word];

		long double valueForWord = 1.0;

		for(ui j = 0; j < wordCount; j++)
		{
			valueForWord *= (baseterm + j);
		}

		finalThirdValueNume *= valueForWord;
		i++;
	}

	double baseTermDenom = dataIOObj.sumAlpha + dataIOObj.NTWSumWordsVec[topic];

	for(ui i = 0; i < doc.size(); i++)
	{
		finalThirdValueDenom *= (baseTermDenom + i);
	}

	finalThirdValue = finalThirdValueNume / finalThirdValueDenom;

	return finalThirdValue;
}

int  GibbsSampler :: printTopicTopicCount(DataIO &dataIOObj)
{
	int sumAll = 0;
	for(ui i = 0; i < dataIOObj.numTopics; i++)
	{
		for(ui j = 0; j < dataIOObj.numTopics; j++)
		{
			// cout << NTTCountVec[i][j] << " ";
			sumAll += dataIOObj.NTTCountVec[i][j];
		}
		// cout << "\n";
	}
	cout << "Sum Over all i-j -- " << sumAll << "\n";
	
	// if(sumAll > countInteractions)
	// {
	// 	cout << "Sum all is greater than Count Interactions (For first iteration it can be...)\n";
	// 	cout << sumAll << " " << countInteractions << "\n";
	// 	// exit(0);
	// }

	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////
///////  SAMPLE PARENT   //////
///////////////////////////////

int GibbsSampler :: sampleParentAssignment(DataIO &dataIOObj, int ITE, int BURN_IN)
{
	cout << "Sampling Parent Assignment ... \n";

	dataIOObj.nodeNodeCombUpdateInfluence.clear();
	// logLikelihood = 0;
	// not using nodeNodeCombUpdateInfluence for now...
	// this is to keep the edges that have change in the parent...

	li previousParentNode, assignedParentNode;

	for(unsigned int i = 0; i < dataIOObj.newSyntheticEvents.size(); i++)
	{
		li eventIndex = dataIOObj.newSyntheticEvents[i][0];
		li eventNode = dataIOObj.newSyntheticEvents[i][1];

		li eventParent = dataIOObj.newSyntheticEvents[i][2];
		li eventTopic = dataIOObj.newSyntheticEvents[i][3];


		double eventTime = dataIOObj.eventIndexTimestamps[i];

		vector<ui> doc;
		// cout << "eventIndex -- " << eventIndex << endl;

		li previousParent = eventParent;

		if(previousParent >= 0)
		{
			previousParentNode = dataIOObj.newSyntheticEvents[previousParent][1];
		}
		else
		{
			previousParentNode = eventNode;
		}

		decreamentCountFromMatrices(dataIOObj, eventIndex, eventNode, eventParent, eventTopic, doc, 0);
		// cout << "Done Decrementing Counters...\n";
		li assignedParent = getSampledParentAssignment(dataIOObj, eventTime, eventNode, eventIndex, eventTopic, ITE, BURN_IN);
		// cout << "Got sampled Parent assignment...-- " << eventIndex << endl;
		increamentCountToMatrices(dataIOObj, eventIndex, eventNode, assignedParent, eventTopic, doc, 0);
		// cout << "Done INcrementing Counters...\n";
		// assign the new sampled parent event to the event...
		
		dataIOObj.newSyntheticEvents[eventIndex][2] = assignedParent;

		if(assignedParent >= 0)
		{
			assignedParentNode = dataIOObj.newSyntheticEvents[assignedParent][1];
		}
		else
		{
			assignedParentNode = eventNode;
		}

		if(previousParent == -1 && assignedParent >= 0)
		{
			dataIOObj.nodeNodeCombUpdateInfluence[assignedParentNode].push_back(eventNode);
		}
		else if(previousParent >= 0 && assignedParent == -1)
		{
			dataIOObj.nodeNodeCombUpdateInfluence[previousParentNode].push_back(eventNode);
		}
		else if(previousParent >= 0 && assignedParent >= 0)
		{
			dataIOObj.nodeNodeCombUpdateInfluence[assignedParentNode].push_back(eventNode);
			dataIOObj.nodeNodeCombUpdateInfluence[previousParentNode].push_back(eventNode);
		}

		// if(i % 10000 == 0)
			// printf("Parent Assignment -- %d\n", i);

		if(ITE > BURN_IN && ITE%10 == 0)
		{
			dataIOObj.parentDistEvery10thIter[eventIndex].push_back(assignedParent);
		}
	}

	// cout << "\nLog Likelihood = " << logLikelihood << "\n";

	return 0;
}

li GibbsSampler :: getSampledParentAssignment(DataIO &dataIOObj, double eventTime, li eventNode, li eventIndex, li eventTopic, int ITE, int BURN_IN)
{
	li assignedParent = -2;
	// cout << "calculatedProbVec-1 -- " << dataIOObj.allPossibleParentEvents[eventIndex].size() << endl;
	vector <ui> possibleParentEvents = dataIOObj.allPossibleParentEvents[eventIndex];
	// cout << "calculatedProbVec-2 -- "<< possibleParentEvents.size()  << "\n";
	
	vector <double> calculatedProbVec;

	if(possibleParentEvents.size() > 0)
	{
		vector<double> possibleParentExp = dataIOObj.allPossibleParentEventsExponentials[eventIndex];
		// cout << "calculatedProbVec-2-exp -- "<< possibleParentExp.size()  << "\n";

		calculatedProbVec = populateCalculatedProbVec(dataIOObj, possibleParentEvents, possibleParentExp, eventNode, eventTopic, eventTime, ITE);
		// cout << "calculatedProbVec-3 \n";
		// cout << "EventIndex -- " << eventIndex << "-- calculatedProbVec.size -- " << calculatedProbVec.size() << endl;

		// ui sampledIndex = getSampleFromMultinomial(calculatedProbVec);
		ui sampledIndex = getSampleFromDiscreteDist(calculatedProbVec);
		// cout << "getSampleFromDiscreteDist \n";

		if(sampledIndex == calculatedProbVec.size()-1)
		{
			assignedParent = -1;
		}
		else
		{
			assignedParent = possibleParentEvents[sampledIndex];
		}

		logLikelihood += log(calculatedProbVec[sampledIndex]);
	}
	else
	{
		assignedParent = -1;
		calculatedProbVec.push_back(1);
	}

	// updating for avgProbability of the parent events... at the end lets write it to a file...  
	if(ITE > BURN_IN && ITE%10 == 0)
	{
		for(ui i = 0; i < dataIOObj.avgProbParForAllEvents[eventIndex].size(); i++)
		{
			dataIOObj.avgProbParForAllEvents[eventIndex][i] += calculatedProbVec[i];
		}
	}

	if(assignedParent > (int)dataIOObj.allEvents.size()-1)
	{
		cout << "some issue with parent sampling...\n";
		exit(0);
	}

	return assignedParent;
}

vector<double> GibbsSampler :: populateCalculatedProbVec(DataIO &dataIOObj, vector <ui> possibleParentEvents, vector<double> possibleParentExp, li eventNode, li eventTopic, double eventTime, int ITE)
{
	// cout << "calculatedProbVec \n";
	// the last index is for probability of self event being the parent...
	vector <double> calculatedProbVec(possibleParentEvents.size() + 1, 0.0);
	
	// double firstTermNume, firstTermDenom, firstTerm, secondTerm;
	double firstTerm, secondTerm;

	double normalizationFactor = 0;

	unsigned int k;
	for(k = 0; k < possibleParentEvents.size(); k++)
	{
		vector <li> possParentEvent = dataIOObj.newSyntheticEvents[possibleParentEvents[k]];
		// li possParentIndex = possParentEvent[0];
		li possParentNode = possParentEvent[1];
		li possParentEventTopic = possParentEvent[3];

		firstTerm = getFirstTermOfParentAssignment(dataIOObj, possParentEventTopic, eventTopic);
		// cout << "getFirstTermOfParentAssignment \n";

		int outDeg = dataIOObj.outDegreeMap[possParentNode];
		int inDeg = dataIOObj.inDegreeMap[eventNode];
		int outDegGID = floor(log2(outDeg)/log2(logBase));
		int inDegGID = floor(log2(inDeg)/log2(logBase));

		string gid = to_string(outDegGID).append("_").append(to_string(inDegGID));

		double Wuv;
		
		// Wuv = userUserInfluence[possParentNode][eventNode];
		Wuv = dataIOObj.groupWuv[gid];

		if(Wuv > 0)
		{
			double temporalDecay = possibleParentExp[k];
			secondTerm = (Wuv * temporalDecay);
		}
		else
		{
			secondTerm = 0;
		}

		// calculatedProbVec[k] = firstTerm + secondTerm;
		calculatedProbVec[k] = firstTerm * secondTerm;

		normalizationFactor += calculatedProbVec[k];
	}

	// get prob of having no parent
	firstTerm = getFirstTermOfParentAssignmentNoParent(dataIOObj, eventNode, eventTopic);
	// cout << "getFirstTermOfParentAssignment \n";
	
	if(dataIOObj.userBaseRateMap.find(eventNode) != dataIOObj.userBaseRateMap.end())
	{
		secondTerm = dataIOObj.userBaseRateMap[eventNode];
	}
	else
	{
		secondTerm = 0;	
		// secondTerm = defaultMuVal;	
	}
	
	// secondTerm = defaultMuVal;	

	calculatedProbVec[k] = firstTerm * secondTerm;

	normalizationFactor += calculatedProbVec[k];

	// it might happen that all the probabilities turn out to be zero... and would have problem with loglikelihood to be nan...
	// Thus making this event as the spontaneous event... and assigning all the probability to self event...
	if(normalizationFactor == 0)
	{
		calculatedProbVec[calculatedProbVec.size() - 1] = 1;
		normalizationFactor = 1;
	}
	
	for(ui i = 0; i < calculatedProbVec.size(); i++)
	{
		calculatedProbVec[i] = calculatedProbVec[i] / normalizationFactor;
	}

	return calculatedProbVec;
}

double GibbsSampler :: getFirstTermOfParentAssignment(DataIO &dataIOObj, li possParentEventTopic, li eventTopic)
{
	double firstTerm = 0.0;
	double firstTermNume = 0.0;
	double firstTermDenom = 1.0;
	
	firstTermNume = dataIOObj.hyperBeta[eventTopic] + dataIOObj.NTTCountVec[possParentEventTopic][eventTopic];
	
	firstTermDenom = dataIOObj.sumBeta + dataIOObj.NTTSumTopicsVec[possParentEventTopic];

	firstTerm = firstTermNume / firstTermDenom;

	return firstTerm;
}

double GibbsSampler :: getFirstTermOfParentAssignmentNoParent(DataIO &dataIOObj, li eventNode, li eventTopic)
{
	// cout << "in getFirstTermOfParentAssignmentNoParent \n";
	double firstTerm = 0.0;
	double firstTermNume = 0.0;
	double firstTermDenom = 1.0;

	firstTermNume = dataIOObj.hyperGamma[eventTopic] + dataIOObj.NUTCountVec[eventNode][eventTopic];
	firstTermDenom = dataIOObj.sumGamma + dataIOObj.NUTSumTopicsVec[eventNode];
	firstTerm = firstTermNume / firstTermDenom;

	return firstTerm;
}

vector<ui> GibbsSampler :: getPossibleParentEvents(DataIO &dataIOObj, li eventNode, li eventIndex)
{
	vector <ui> possibleParentEvents = dataIOObj.allPossibleParentEvents[eventIndex];  
	return possibleParentEvents;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// Decreament Counters
int GibbsSampler :: decreamentCountFromMatrices(DataIO &dataIOObj, li eventIndex, li eventNode, li eventParent, li eventTopic, vector <ui> doc, bool topicSampling)
{
	// cout << "decreament count from topic-topic matrix.. from the cell parentTopic -> eventTopic";
	// decreament only if the eventTopic is valid...
	if(eventTopic >= 0)
	{
		int eventParentTopic;

		if(eventParent >= 0)
		{
			eventParentTopic = dataIOObj.newSyntheticEvents[eventParent][3];
			// printUnorderedMap("ttopic");
			decreamentCountFromTopicTopic(dataIOObj, eventParentTopic, eventTopic);
			// printUnorderedMap("ttopic");
		}
		else if(eventParent == -1)
		{
			// printUnorderedMap("utopic");
			decreamentCountFromUserTopic(dataIOObj, eventNode, eventTopic);
			// printUnorderedMap("utopic");
		}

		if(topicSampling == 1)
		{
			// decreament the count of words corresponding to this event...
			decreamentCountFromTopicWord(dataIOObj, eventTopic, doc);
			
			// decreament counts corresponding to child events for topic sampling...
			decreamentCountsFromChildEvents(dataIOObj, eventTopic, eventIndex);
		}
	}

	return 0;
}


int GibbsSampler :: decreamentCountFromTopicTopic(DataIO &dataIOObj, li eventParentTopic, li eventTopic)
{
	// cout << "decreament count from topic-topic matrix.. from the cell parentTopic -> eventTopic";

	if(dataIOObj.NTTCountVec[eventParentTopic][eventTopic] > 0)
	{
		dataIOObj.NTTCountVec[eventParentTopic][eventTopic] -= 1;
	}

	if(dataIOObj.NTTSumTopicsVec[eventParentTopic] > 0)
	{
		dataIOObj.NTTSumTopicsVec[eventParentTopic] -= 1;
	}

	return 0;
}


int GibbsSampler :: decreamentCountFromUserTopic(DataIO &dataIOObj, li eventNode, li eventTopic)
{
	// cout << "decreament count from user-topic matrix.. from the cell user -> eventTopic";

	if(dataIOObj.NUTCountVec[eventNode][eventTopic] > 0)
	{
		dataIOObj.NUTCountVec[eventNode][eventTopic] -= 1;	
	}

	if(dataIOObj.NUTSumTopicsVec[eventNode] > 0)
	{
		dataIOObj.NUTSumTopicsVec[eventNode] -= 1;	
	}

	return 0;
}

int GibbsSampler :: decreamentCountFromTopicWord(DataIO &dataIOObj, li eventTopic, vector <ui> doc)
{
	if(dataIOObj.NTWSumWordsVec[eventTopic] >= doc.size())
	{
		dataIOObj.NTWSumWordsVec[eventTopic] -= doc.size();
	}
	else
	{
		// cout << "Somethings wrong with the Topic Word Counts... \n";
		dataIOObj.NTWSumWordsVec[eventTopic] = 0;
	}

	for(unsigned int i = 0; i < doc.size(); i++)
	{
		if(dataIOObj.NTWCountVec[eventTopic][doc[i]] > 0)
		{
			dataIOObj.NTWCountVec[eventTopic][doc[i]] -= 1;
		}
	}

	return 0;
}


int GibbsSampler :: decreamentCountsFromChildEvents(DataIO &dataIOObj, li eventTopic, li eventIndex)
{
	// vector <li> childEventsList;
	vector <ui> childEventsList;
	try
	{
		childEventsList =  dataIOObj.childEventsMap.at(eventIndex);
	}
	catch(exception &e)
	{
		return 0;
	}

	if(childEventsList.size() > 0)
	{
		for(unsigned int i = 0; i < childEventsList.size(); i++)
		{
			li childEventIndex = childEventsList[i];

			li childEventTopic = dataIOObj.newSyntheticEvents[childEventIndex][3];

			if(childEventTopic != -1)
			{
				decreamentCountFromTopicTopic(dataIOObj, eventTopic, childEventTopic);
			}
		}
	}

	return 0;
}


// Increament Counters
int GibbsSampler :: increamentCountToMatrices(DataIO &dataIOObj, li eventIndex, li eventNode, li eventParent, li eventTopic, vector<ui> doc, bool topicSampling)
{
	// cout << "increament count to topic-topic matrix... in the cell parentTopic -> assignedTopic";

	if(eventParent >= 0)
	{
		int eventParentTopic = dataIOObj.newSyntheticEvents[eventParent][3];
		// printUnorderedMap("ttopic");
		// cout << "increamenting...\n";
		increamentCountInTopicTopic(dataIOObj, eventParentTopic, eventTopic);
		// printUnorderedMap("ttopic");
	}
	else if (eventParent == -1)
	{
		// printUnorderedMap("utopic");
		// cout << "increamenting...\n";
		increamentCountInUserTopic(dataIOObj, eventNode, eventTopic);
		// printUnorderedMap("utopic");
	}

	if(topicSampling == 1)
	{
		// increament counts in topic-word...
		increamentCountInTopicWord(dataIOObj, eventTopic, doc);
		// printUnorderedMap("psitopic");
		
		// increament counts corresponding to child events...
		increamentCountsForChildEvents(dataIOObj, eventTopic, eventIndex);
	}

	return 0;
}

int GibbsSampler :: increamentCountInTopicTopic(DataIO &dataIOObj, li eventParentTopic, li eventTopic)
{

	dataIOObj.NTTSumTopicsVec[eventParentTopic] += 1;
	dataIOObj.NTTCountVec[eventParentTopic][eventTopic] += 1;

	if(dataIOObj.NTTCountVec[eventParentTopic][eventTopic] > dataIOObj.newSyntheticEvents.size())
	{
		cout << "Somethings wrong.. the size of NTT is more than number of events...\n";
		cout << dataIOObj.NTTCountVec[eventParentTopic][eventTopic] << " " << dataIOObj.newSyntheticEvents.size() << "\n";
		exit(0);
	}


	if(dataIOObj.NTTSumTopicsVec[eventParentTopic] > dataIOObj.newSyntheticEvents.size())
	{
		cout << "Somethings wrong.. the size of NTT (sum) is more than number of events...\n";
		cout << dataIOObj.NTTSumTopicsVec[eventParentTopic] << " " << dataIOObj.newSyntheticEvents.size() << "\n";
		exit(0);
	}

	// NTSumTopics[eventParentTopic] += 1;

	// li cellIndex = getCellKey(eventParentTopic, eventTopic);
	// NTTCountMatrix[cellIndex] += 1;

	return 0;
}

int GibbsSampler :: increamentCountInUserTopic(DataIO &dataIOObj, li eventNode, li eventTopic)
{
	dataIOObj.NUTSumTopicsVec[eventNode] += 1;
	dataIOObj.NUTCountVec[eventNode][eventTopic] += 1;
	
	if(dataIOObj.NUTCountVec[eventNode][eventTopic] > dataIOObj.newSyntheticEvents.size())
	{
		cout << "Somethings wrong.. the size of NUT is more than number of events...\n";
		cout << dataIOObj.NUTCountVec[eventNode][eventTopic] << " " << dataIOObj.newSyntheticEvents.size() << "\n";
		exit(0);
	}


	if(dataIOObj.NUTSumTopicsVec[eventNode] > dataIOObj.newSyntheticEvents.size())
	{
		cout << "Somethings wrong.. the size of NUT (sum) is more than number of events...\n";
		cout << dataIOObj.NUTSumTopicsVec[eventNode] << " " << dataIOObj.newSyntheticEvents.size() << "\n";
		exit(0);
	}


	// NUSumTopics[eventNode] += 1;

	// li cellIndex = getCellKey(eventNode, eventTopic);
	// NUTCountMatrix[cellIndex] += 1;

	return 0;
}

int GibbsSampler :: increamentCountInTopicWord(DataIO &dataIOObj, li eventTopic, vector<ui> doc)
{
	dataIOObj.NTWSumWordsVec[eventTopic] += doc.size();
	// NTSumWords[eventTopic] += doc.size();

	if(dataIOObj.NTWSumWordsVec[eventTopic] > dataIOObj.totalWords)
	{
		cout << "Total words = " << dataIOObj.totalWords << " NTW sum = " << dataIOObj.NTWSumWordsVec[eventTopic] << "\n";
		exit(0);
	}

	for(unsigned int i = 0; i < doc.size(); i++)
	{
		dataIOObj.NTWCountVec[eventTopic][doc[i]] += 1;

		if(dataIOObj.NTWCountVec[eventTopic][doc[i]] > dataIOObj.totalWords)
		{
			cout << "Total words = " << dataIOObj.totalWords << " NTW = " << dataIOObj.NTWCountVec[eventTopic][doc[i]] << " " << doc[i] << "\n";
			exit(0);
		}
	}
	
	return 0;
}

int GibbsSampler :: increamentCountsForChildEvents(DataIO &dataIOObj, li eventTopic, li eventIndex)
{
	// vector <li> childEventsList;
	vector <ui> childEventsList;
	try
	{
		childEventsList = dataIOObj.childEventsMap.at(eventIndex);
	}
	catch(exception &e)
	{
		return 0;
	}

	if(childEventsList.size() > 0)
	{
		// if(childEventsList.size() > 42)
			// printf("ChildEvents.size() -- %lu -- EventIndex -- %lu\n", childEventsList.size(), eventIndex);

		for(unsigned int i = 0; i < childEventsList.size(); i++)
		{
			li childEventIndex = childEventsList[i];

			li childEventTopic = dataIOObj.newSyntheticEvents[childEventIndex][3];

			if(childEventTopic != -1)
			{
				increamentCountInTopicTopic(dataIOObj, eventTopic, childEventTopic);
			}
		}
	}

	return 0;
}


int GibbsSampler :: getSampleFromMultinomial(vector<double> calculatedProbVec)
{
	int assignedInd;

	// get normalized prob vector...
	vector <double> normalizedProbVector = getNormalizedLogProb(calculatedProbVec);
	
	assignedInd = getSampleFromDiscreteDist(normalizedProbVector);

	logLikelihood += log(normalizedProbVector[assignedInd]);

    return assignedInd;
}