#include "namespace.h"
#include "dataIO.h"
#include "initialization.h"
#include "utilities.h"

using namespace std;

InitializeModel :: InitializeModel(DataIO &dataIOObj)
{
	baseAlpha = 0.01;
	baseBeta = 1;
	logBase = 1.5;

	initializeForSampler(dataIOObj);

	updateNodeNodeCountMap(dataIOObj);
	
	initializeUserUserInfluence(dataIOObj);

	initializeAvgProbabilityVectors(dataIOObj);
	initializeAvgTopicProbabilityVectors(dataIOObj);

	cout << "Done with Initializing...\n";
}


// vector< vector <li> > InitializeModel :: initializeForSampler(vector< vector <li> > allEvents, int flag)
int InitializeModel :: initializeForSampler(DataIO &dataIOObj)
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

int InitializeModel :: updateCountMatrices(DataIO &dataIOObj, vector< vector <li> > &localNewSyntheticEvents)
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


int InitializeModel :: updateNodeNodeCountMap(DataIO &dataIOObj)
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



int InitializeModel :: initializeUserUserInfluence(DataIO &dataIOObj)
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



int InitializeModel ::  initializeAvgProbabilityVectors(DataIO &dataIOObj)
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

int InitializeModel :: initializeAvgTopicProbabilityVectors(DataIO &dataIOObj)
{
	for(ui i = 0; i < dataIOObj.allDocs.size(); i++)
	{
		vector <double> initVec(dataIOObj.numTopics, 0.0);
		dataIOObj.avgTopicProbVector.push_back(initVec);
		initVec.clear();
	}
	return 0;
}