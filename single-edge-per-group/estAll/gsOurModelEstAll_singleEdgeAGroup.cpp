#include <iostream>				//for basic C++ functions 
#include <cstdio>				//for std C functions
#include <fstream>				//for i/o stream
#include <sstream>				//for parse data using stringstreams...
#include <string>				//for C++ string functions
#include <cstring>				//we need this for memset... and string functions from C
#include <unordered_map>		//to maintain the map of user-tweet_count, and user_user_tweet_count
#include <map>		//to maintain the map of user-tweet_count, and user_user_tweet_count
#include <ctime>				//for time realted functions... to convert user-formatted time to timestamp
#include <vector>
#include <assert.h>				//for assertions
#include <limits>		
#include <cmath>
// #include <math.h>
// #include <cdouble>
#include <algorithm>
#include <numeric>	
#include <random>
#include <chrono>
#include <cfenv>
#include <iomanip>


// #pragma STDC FENV_ACCESS ON

using namespace std;
using namespace std::chrono;


// using ull = unsigned long long int;
using ll = long long int;
using li = long int;
using ui = unsigned int;

#define INF std::numeric_limits<int>::max()					//use limits for infinity
#define PI 3.14159265

// function declaration

int getTopicTopicCombinations();
int printTopicTopicCount();

// Read data from files....
// vector < vector <double> > getSyntheticEventsFromFile(string fileName);
vector < vector <li> > getSyntheticEventsFromFile(string fileName);
vector < vector <ui> > getSyntheticDocsFromFile(string fileName);
// unordered_map <ll, vector <ll> > readIntVectorMapFromFile(string fileName);
vector < vector <ui> > readIntVectorMapFromFile(string fileName);
vector< vector <double> > readMultipleDoubleVectorsFromFile(string fileName);
vector <double>  readDoubleVectorFromFile(string fileName);
unordered_map <li, unordered_map <li, double> > readUserUserInfluenceFromFile(string fileName);
unordered_map <ui, unordered_map <ui, float> > getUserUserInfluence(string fileName);
map<ui, double> getUserBaseRates(string fileName);
int readHyperParametersVectorsFromFile();

// initialization
vector< vector <li> > initializeForSampler(vector< vector <li> > allEvents, int flag);
int initializeUserUserInfluence();
double getPriorRateParamFromData();
int initializeAndUpdateAlphaValue();
int updateUserBaseRates();
int initializeAvgProbabilityVectors();
int initializeAvgTopicProbabilityVectors();

// sample topic
int sampleTopicAssignment(int ITE);
li getSampledTopicAssignment(li eventIndex, li eventNode, li eventParent, vector<ui> doc, int ITE);
double getFirstTermOfTopicAssignmentCondProb(li eventNode, li eventParent, li topic);
double getMiddleTermOfTopicAssignmentCondProb(li topic, unordered_map <int, ui> childEventTopicsHist, li eventIndex, li eventParent);
double getThirdTermOfTopicAssignmentCondProb(li topic, vector<ui> doc, vector<ui> wordHistVec);
unordered_map <int, ui> getHistOfTopicsOverChildEvents(li eventIndex);
unordered_map <ui, ui> getHistOfWordsOverWordsFromDoc(vector<ui> doc);
int createWordHistForAllDocs();


// sample parent
int sampleParentAssignment(int ITE);
li getSampledParentAssignment(double eventTime, li eventNode, li eventIndex, li eventTopic, int ITE);
vector<double> populateCalculatedProbVec(vector <ui> possibleParentEvents, vector<double> possibleParentExp, li eventNode, li eventTopic, double eventTime, int ITE);
double getFirstTermOfParentAssignment(li possParentEventTopic, li eventTopic);
double getFirstTermOfParentAssignmentNoParent(li eventNode, li eventTopic);
vector<ui> getPossibleParentEvents(li eventNode, li eventIndex);
vector<double> computeParentExponentials(vector<ui> possibleParentEvents, li eventIndex);


// sample Influence
int sampleInfluenceAssignment(int ITE);
int updateNodeNodeCountMap();


// decreament
int decreamentCountFromMatrices(li eventIndex, li eventNode, li eventParent, li eventTopic, vector <ui> doc, bool topicSampling);
int decreamentCountFromTopicTopic(li eventParentTopic, li eventTopic);
int decreamentCountFromUserTopic(li eventNode, li eventTopic);
int decreamentCountFromTopicWord(li eventTopic, vector <ui> doc);
int decreamentCountsFromChildEvents(li eventTopic, li eventIndex);


// increament
int increamentCountToMatrices(li eventIndex, li eventNode, li eventParent, li eventTopic, vector<ui> doc, bool topicSampling);
int increamentCountInTopicTopic(li eventParentTopic, li eventTopic);
int increamentCountInUserTopic(li eventNode, li eventTopic);
int increamentCountInTopicWord(li eventTopic, vector<ui> doc);
int increamentCountsForChildEvents(li eventTopic, li eventIndex);

// validations
int countCorrectFractionTopicAssignments();
int countCorrectFractionParentAssignments();
int writeEvery10ItersParentAssignmentToFile();
int writeEvery10ItersTopicAssignmentToFile();
int writeAvgProbVectorsToFile();
int writeAvgTopicProbVectorsToFile();
int writeTopicTopicInteractionToFile();
int writeUserUserInfluenceToFile();
int writeUserBaseRatesToFile();
int writeNodeNodeCountAndNodeCount();

// util functions
li getCellKey(li firstPart, li secondPart);
int getSampleFromMultinomial(vector<double> calculatedProbVec);
vector<double> getNormalizedLogProb(vector<double> calculatedProbVec);
int getSampleFromDiscreteDist(vector<double> normalizedProbVector);
double getSampleFromGamma(double alpha, double beta);
int printUnorderedMap(string matrixName);
string& SSS (const char* s);


int getAvgChildEvents();
int populateParentEventsForAll();
int populateParentEventsForAllFromFile(string parentsFile, string parentExpFile);
int writeAssignmentsToFile();

int printTopicTopicCorrelation();

bool sortcol(const vector<float>& v1, const vector<float>& v2)
{
	// return v1[1] < v2[1];
	return v1[1] > v2[1]; 				// sort descending...
}


// global members declaration
vector< vector <li> > allEvents;
vector< vector <li> > newSyntheticEvents;
vector< vector <ui> > allDocs;

// unordered_map <li, unordered_map<li, double> > origUserUserInfluence;
// unordered_map <li, unordered_map<li, double> > userUserInfluence;

vector<double> hyperBeta, hyperGamma, hyperAlpha;
double sumAlpha, sumBeta, sumGamma;
vector< vector <double> > topicTopicProbVector,  userTopicPrefVector, topicWordProbVector;

// to maintain a map of eventIndex and docSize...
unordered_map <ui, int> eventIndexDocSizeMap;

map<ui, vector<ui> > nodeEventsMap;

// alpha, beta gamma initialization for userUserInfluence...
double baseAlpha = 0.01, baseBeta = 1;
double rateParam;

// double muVal = 0.00002;
double defaultMuVal = 2.16428e-06;

map<ui, double> userBaseRateMap;

// double timeKernelMultiplier = 0.01;
double timeKernelMultiplier = 1;

ui maxNumNodes = 1000002;					// 7697889				
// ui maxDegree = 600;
// unsigned maxDegree = 248438;

int BURN_IN = 200;
ui ITERATIONS = 301;

ui numNodes = maxNumNodes;

ui numTopics = 100;

// ui vocabSize = 64500; 						// 10000
ui vocabSize = 33900; 							// 10000
ui totalWords = 0;

double alphaWord = 0.1;
// double sumAlpha;

// Count Matrices... 
// std::vector<std::vector<int>> vec_2d(rows, std::vector<int>(cols, 0));
vector < vector < ui > > NTWCountVec(numTopics, vector < ui >(vocabSize, 0));
vector < vector < ui > > NTTCountVec(numTopics, vector < ui >(numTopics, 0));
vector < vector < ui > > NUTCountVec(maxNumNodes, vector < ui >(numTopics, 0));

vector < ui > NTWSumWordsVec(numTopics);
vector < ui > NTTSumTopicsVec(numTopics);
vector < ui > NUTSumTopicsVec(numNodes);

vector < vector <ui> > wordHistAllDocsVector;

// unordered_map <li, unordered_map<li, double> > origUserUserInfluence;
unordered_map <ui, unordered_map<ui, double> > userUserInfluence;
// vector < vector <ll> > followersMap(maxNumNodes, vector < ll >(maxDegree + 1, -1) );
vector < vector <ui> > followersMap;
// vector < vector <ll> > reverseFollowersMap(maxNumNodes, vector < ll >(maxDegree + 1, -1) );
vector < vector <ui> > reverseFollowersMap;

vector <vector <ui> > allPossibleParentEvents;
vector <vector <double> > allPossibleParentEventsExponentials;
vector <vector <double> > avgProbParForAllEvents;
vector <vector <double> > avgTopicProbVector;


// unordered_map<ui, ui> nodeEventsCountMap;
vector<ui> nodeEventsCountMap(maxNumNodes, 0);
unordered_map<ui, unordered_map<ui, ui> > nodeNodeCount;
unordered_map<ui, double> eventIndexTimestamps;

// unordered_map <ui, vector<ui> > childEventsMap;
map <ui, vector<ui> > childEventsMap;
vector<ui> levelInfo;

unordered_map<ui, vector <ui> > nodeNodeCombUpdateInfluence;

// validation Maps
// unordered_map <int, vector<int> > topicDistEvery10thIter, parentDistEvery10thIter;
map <int, vector<int> > topicDistEvery10thIter, parentDistEvery10thIter;

int flag = 0;
int ofHitCount = 0;

ui maxLevel;

long int countInteractions = 0;
long int totalDecreaments = 0;
long int totalIncreaments = 0;

double logLikelihood = 0;
double alphaEstimate = 0;

int main(int argc, char *argv[])
{

	if(argc > 1)
	{
		BURN_IN = atoi(argv[1]);
		ITERATIONS = atoi(argv[2]);
	}
	
	cout << "Will be running sampler for BURN_IN = " << BURN_IN << " And ITERATIONS = " << ITERATIONS << "\n";

	// allEvents = getSyntheticEventsFromFile("../centralFiles/allSyntheticEvents_Orig.txt");
	// allEvents = getSyntheticEventsFromFile("./allSemiSyntheticEvents_global_C5_AlphaExp_1.0.txt");
	cout << "Reading Events...\n";
	allEvents = getSyntheticEventsFromFile("../events_semisyn_sample_gamma_1M.txt");
	// allEvents = getSyntheticEventsFromFile("../events_semisyn_sample.txt");
	cout << "Got all the events... " << allEvents.size() << "\n";

	allDocs = getSyntheticDocsFromFile("../docs_semisyn_sample_gamma_1M.txt");
	// allDocs = getSyntheticDocsFromFile("../docs_semisyn_sample.txt");
	cout << "Got all the docs...\n";

	cout << "Creating hist of words for each doc...\n";
	createWordHistForAllDocs();
	cout << "Created hist of words for each doc...\n";

	cout << "Getting followers map\n";
	followersMap = readIntVectorMapFromFile("../../mapped_users_followers_restricted_tweets_top_5k_hashes.txt");			// this would be required for the Wuv matrix.... 
	cout << "Got the followers map... " << followersMap.size() << "\n";

	readHyperParametersVectorsFromFile();	// for now not required...

	// populate all possible parent events...
	cout << "Getting possible parent events for all the events...\n";
	populateParentEventsForAllFromFile("../top100CandParents_sample_events_1L_time_15.txt", "../top100CandParExp_sample_events_1L_time_15.txt");
	// populateParentEventsForAllFromFile("../top100CandParents_sample_events_allevents.txt", "../top100CandParExp_sample_events_allevents.txt");

	cout << "Got possible parent events for all the events..." << allPossibleParentEvents.size() << " " << allPossibleParentEventsExponentials.size() << "\n";

	newSyntheticEvents = initializeForSampler(allEvents, flag);
	// newSyntheticEvents = allEvents;

	// Not using this...
	// rateParam = getPriorRateParamFromData();
	// cout << "Rate param = " << rateParam << "\n";

	// initialize the influence matrix to the or some small value...
	initializeUserUserInfluence();
	// As we have the initial (time based) parent assignment, we can get some estimate of the Wuv values...
	updateNodeNodeCountMap();
	// We can get some initial Wuv estimate...
	sampleInfluenceAssignment(-1);

	// update the user base rates... we can get some estimate of the base rates from here...
	updateUserBaseRates();
	
	// cout << "Getting user base rates...\n";
	// userBaseRateMap = getUserBaseRates("./oscars_users_Baserate_previousNeighborParent.txt");
	// userBaseRateMap = getUserBaseRates("./oscars_users_baserate.txt");							// we might need this while assigning self parent.. or we take some default value...

	// cout << "Total Number of interactions = " << countInteractions << "\n";


	initializeAvgProbabilityVectors();
	initializeAvgTopicProbabilityVectors();				

	cout << "done with initialization ... \n";

	ofstream llFile;
	llFile.open("estAllLogLikelihood.txt");
	
	
	for(ui ITE = 0; ITE < ITERATIONS; ITE++)
	{
		logLikelihood = 0;

		high_resolution_clock::time_point t1, t2;
		auto duration = 0;

		// we are not dealing with topic assignment... In this case the topic assignments are known...
		 
		t1 = high_resolution_clock::now();
		sampleTopicAssignment(ITE);
		t2 = high_resolution_clock::now();
	
		duration = duration_cast<microseconds>( t2 - t1 ).count();
		cout << "sampled topics for iteration -- "<< ITE << "--  Time taken -- " << duration << endl;


		t1 = high_resolution_clock::now();
		sampleParentAssignment(ITE);
		t2 = high_resolution_clock::now();
	
		duration = duration_cast<microseconds>( t2 - t1 ).count();
		cout << "sampled parents for iteration -- "<< ITE << "--  Time taken -- " << duration << endl;
	

		// updating the childEventsMap... as it might change after each iteration once we incorporate parent sampling
		t1 = high_resolution_clock::now();
		// As there is the change in the parent events, update the child events...
		childEventsMap.clear();
		for(unsigned int i = 0; i < newSyntheticEvents.size(); i++)
		{
			li eventCurrParent = newSyntheticEvents[i][2];
			if(eventCurrParent != -1)
			{
				childEventsMap[eventCurrParent].push_back(i);
			}
		}
		// As there is the change in the parent assignments... node-node counts (Nuv) would change over edges...
		updateNodeNodeCountMap();
	
		t2 = high_resolution_clock::now();
		duration = duration_cast<microseconds>( t2 - t1 ).count();
		cout << "Updating Child Events Map and Node-Node Count after -- "<< ITE << " Iteration --  Time  taken -- " << duration << endl;

		// calculate the fresh Wuv values with the update Node-Node counts...
		t1 = high_resolution_clock::now();
		sampleInfluenceAssignment(ITE);
		t2 = high_resolution_clock::now();
	
		duration = duration_cast<microseconds>( t2 - t1 ).count();
		cout << "Sampled User-User Influence for iteration -- "<< ITE << "--  Time taken -- " << duration << endl;

		t1 = high_resolution_clock::now();
		updateUserBaseRates();
		t2 = high_resolution_clock::now();
	
		duration = duration_cast<microseconds>( t2 - t1 ).count();
		cout << "Updating Alpha Value and User Base Rates -- "<< ITE << "--  Time taken -- " << duration << endl;
        cout << "Default Mu Val = " << defaultMuVal << endl;

		// cout << "Log Likelihood = " << logLikelihood << "\n";
		// llFile << logLikelihood << "\n";

		cout << setprecision(10) << "Log Likelihood = " << logLikelihood << "\n";
		llFile << setprecision(10) << logLikelihood << "\n";

		cout << "ITE -- " << ITE << endl;
		
		if(ITE % 100 == 0)
		{
			cout << "ITE -- " << ITE << endl;
		}
	}	// ITERATIONS...

	cout << "Done... Will be writing results to files...\n";

	// cout << "OverFlow Hit Count = " << ofHitCount << endl;
	llFile.close();

	writeEvery10ItersTopicAssignmentToFile();
	writeEvery10ItersParentAssignmentToFile();
	writeAvgProbVectorsToFile();
	writeAvgTopicProbVectorsToFile();
	writeTopicTopicInteractionToFile();				// we will build this from avg probability parent assignment...
	writeUserUserInfluenceToFile();
	writeUserBaseRatesToFile();
	writeNodeNodeCountAndNodeCount();

	return 0;
}

///////////////////////////////
///////// INITIALIZATION //////
///////////////////////////////

vector< vector <li> > initializeForSampler(vector< vector <li> > allEvents, int flag)
{
	// Random topic initialization...
	// Parent Initialization is time based parent...
	vector< vector <li> > localNewSyntheticEvents;
	
	for(ui i = 0; i < allEvents.size(); i++)
	{
		vector<li> tempEvent;
		tempEvent.push_back(i);
		tempEvent.push_back(allEvents[i][1]);

		int currentEventNode = allEvents[i][1];

		if(allPossibleParentEvents[i].size() > 0)
		{
			// initialize to nearest in time parent...
			int possParentEventId = allPossibleParentEvents[i][0];

			tempEvent.push_back(possParentEventId);				// some initialization at the start.... i.e. nearest in time from neighbor...

			// populate the childevents for the parent events...
			childEventsMap[possParentEventId].push_back(i);

			int possParentEventNode = localNewSyntheticEvents[possParentEventId][1];
			// checking or updating the map that keeps track of the updated edges...
			nodeNodeCombUpdateInfluence[possParentEventNode].push_back(currentEventNode);
			// getting the Nuv value... will be updated after every iteration...
			nodeNodeCount[possParentEventNode][currentEventNode]++;
		}
		else
		{
			// if no nearest in time parent, assign no parent...
			tempEvent.push_back(-1);
		}
		// assigning random topic...
		tempEvent.push_back(i % numTopics);										// for some arbit initialization...

		localNewSyntheticEvents.push_back(tempEvent);
		
		// populate the node and eventIndex map... will be useful while searching for possible parent events...
		nodeEventsMap[(ui)currentEventNode].push_back(i);
		// the events here must be already in sorted order, because of the way they are added....

		// i = eventIndex...
		nodeEventsCountMap[(ui)currentEventNode]++;
		// maintaining count of events corresponding to each node ... this will not change irrespective of the gibbs updates... 
		// maintaining a map of eventIndex and timestamp...
		// eventIndexTimestamps[i] = allEvents[i][0];

		// eventIndex++;
	}

	// initialization of the NTT Matrix... because of some initialization of parents...
	int sumAll = 0;
	for(ui i = 0; i < localNewSyntheticEvents.size(); i++)
	{
		int eParent = localNewSyntheticEvents[i][2];
		int eTopic = localNewSyntheticEvents[i][3];

		if(eParent > -1)
		{
			int eParentTopic = localNewSyntheticEvents[eParent][3];
			NTTCountVec[eParentTopic][eTopic] += 1;
			NTTSumTopicsVec[eParentTopic] += 1;
			sumAll++;
		}
		else
		{
			int eNode = localNewSyntheticEvents[i][1];
			NUTCountVec[eNode][eTopic] += 1;
			NUTSumTopicsVec[eNode] += 1;
			// this should also give an estimate of the events without parent...
		}
	}

	cout << "Number of Edges = " << sumAll << "\n";

	// populate the topic-word matrix as well as per the topic assignment...
	for(unsigned int i = 0; i < allDocs.size(); i++)
	{
		// assignedSampleTopics.push_back(-1);
		int initTopic = localNewSyntheticEvents[i][3];;

		for(ui j = 0; j < allDocs[i].size(); j++)
		{
			int word = allDocs[i][j];
			NTWCountVec[initTopic][word]++;
		}

		NTWSumWordsVec[initTopic] += allDocs[i].size();
	}

	return localNewSyntheticEvents;
}


double getPriorRateParamFromData()
{
	vector<ui> allNuValues = nodeEventsCountMap;

	// take only non-zero values from nodeEventsCountMap
	sort(allNuValues.begin(), allNuValues.end());

	vector<ui> nonZeroNuValues;

	for(ui i = 0; i < allNuValues.size(); i++)
	{
		if(allNuValues[i] > 0)
		{
			nonZeroNuValues.push_back(allNuValues[i]);
		}
	}
	cout << "Nodes having atleast one event = " << nonZeroNuValues.size() << "\n";

	// int midInd = nonZeroNuValues.size() / 2;

	// double medianVal = nonZeroNuValues[midInd] * 1.0;

	double maxVal = nonZeroNuValues[nonZeroNuValues.size() - 1] * 1.0;

	// return medianVal;
	return maxVal;
}

int initializeUserUserInfluence()
{
	cout << "Initializing user user influence to a prior...\n";

	// unordered_map<ui, vector <ui> >::iterator followersMapIt;
	ui count = 0;
	ui uNode, vNode;
	for(ui i = 0; i < followersMap.size(); i++)
	{
		vector <ui> followers = followersMap[i];

		uNode = i;
		
		unordered_map <ui, double> tempNodeInf;

		for(ui j = 0; j < followers.size(); j++)
		{
			vNode = followers[j];
			double scaleParam = 1/(baseBeta);

			// tempNodeInf[vNode] = getSampleFromGamma(baseAlpha, scaleParam);
			tempNodeInf[vNode] = baseAlpha * scaleParam;

			// tempNodeInf[vNode] = 0.010;
		}

		userUserInfluence[uNode] = tempNodeInf;

		tempNodeInf.clear();

		count++;

		if(count % 100000 == 0)
			cout << "Initialized for " <<  count << " users...\n";
	}

	cout << "Initialized user user influence to a prior...\n";

	return 0;
}

int initializeAvgProbabilityVectors()
{
	for(ui i = 0; i < allPossibleParentEvents.size(); i++)
	{
		int probVecSize = allPossibleParentEvents[i].size();
		vector<double> initVec(probVecSize + 1, 0.0);
		avgProbParForAllEvents.push_back(initVec);
		initVec.clear();
	}

	return 0;
}

int initializeAvgTopicProbabilityVectors()
{
	for(ui i = 0; i < allDocs.size(); i++)
	{
		vector <double> initVec(numTopics, 0.0);
		avgTopicProbVector.push_back(initVec);
		initVec.clear();
	}
	return 0;
}


///////////////////////////////
/////// User Base Rates  //////
///////////////////////////////


int updateUserBaseRates()
{
	int sponCount = 0;
	int minTime = eventIndexTimestamps[0];
	int maxTime = eventIndexTimestamps[eventIndexTimestamps.size() - 1];

	int totalDataObservedTime = maxTime - minTime;

	map<int, int> eachNodeSponCount;


	map<int, bool> distinctNodes;

	for(ui i = 0; i < newSyntheticEvents.size(); i++)
	{
		if(newSyntheticEvents[i][2] == -1)
		{
			sponCount++;
			distinctNodes[newSyntheticEvents[i][1]] = 1;
			eachNodeSponCount[newSyntheticEvents[i][1]]++;
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

		userBaseRateMap[nodeId] = (nodeSponCount * 1.0) / totalDataObservedTime;

		averageMuVal +=  userBaseRateMap[nodeId];

		if(userBaseRateMap[nodeId] > maxMuVal)
		{
			maxMuVal = userBaseRateMap[nodeId];
		}

		if(userBaseRateMap[nodeId] < minMuVal)
		{
			minMuVal = userBaseRateMap[nodeId];
		}
	}

	defaultMuVal = averageMuVal / eachNodeSponCount.size();
	
	cout << "sponCount = " << sponCount << " minMuVal = " << minMuVal << " maxMuVal = " << maxMuVal << " MLE MuVal = " << defaultMuVal << endl;

	return 0;
}


///////////////////////////////
/////// SAMPLE INFLUENCE //////
///////////////////////////////

/*
int sampleInfluenceAssignment(int ITE)
{
	cout << "Sample User User Influence\n";
	ui Nuv, Nu;
	double alphaDash, betaDash;
	ui uNode, vNode;

	unordered_map<ui, vector<ui> >::iterator nodeNodeUpdateIterator;

	double avgVal = 0;
	int count = 0;

	cout << "nodeNodeUpdates before clearance = " << nodeNodeCombUpdateInfluence.size() << "\n";

	for(nodeNodeUpdateIterator = nodeNodeCombUpdateInfluence.begin(); nodeNodeUpdateIterator != nodeNodeCombUpdateInfluence.end(); nodeNodeUpdateIterator++)
	{
		uNode = nodeNodeUpdateIterator->first;
		vector<ui> updateForNodes = nodeNodeUpdateIterator->second;

		Nu = nodeEventsCountMap[uNode];

		for(ui i = 0; i < updateForNodes.size(); i++)
		{
			vNode = updateForNodes[i];
			Nuv = nodeNodeCount[uNode][vNode];

			alphaDash = Nuv + baseAlpha;
			// betaDash = 1/(Nu + baseBeta);
			// double scaleParam = 1/(Nu + rateParam);
			double scaleParam = 1/(Nu + baseBeta);
			// betaDash = 1.0/Nu;                 			//remove

			// userUserInfluence[uNode][vNode] = getSampleFromGamma(alphaDash, scaleParam);
			userUserInfluence[uNode][vNode] = alphaDash * scaleParam;

			avgVal += userUserInfluence[uNode][vNode];
			count++;
		}

		nodeNodeCombUpdateInfluence[uNode].clear();
	}

	cout << "Avg Wuv = " << avgVal / count << " Update edges = " << count <<  endl;
	nodeNodeCombUpdateInfluence.clear();
	cout << "nodeNodeUpdates cleared count = " << nodeNodeCombUpdateInfluence.size() << "\n";
	return 0;
}
*/

int sampleInfluenceAssignment(int ITE)
{
	unordered_map<ui, unordered_map<ui, ui> >::iterator nodeNodeCountIterator;
	double avgVal = 0;
	int count = 0;

	// get the Nuv counts and calculate the Wuv values... using the prior alphaDash and betaDash 
	for(nodeNodeCountIterator = nodeNodeCount.begin(); nodeNodeCountIterator != nodeNodeCount.end(); nodeNodeCountIterator++)
	{
		ui uNode = nodeNodeCountIterator->first;
		// Node events count map does not change over iterations...
		ui Nu = nodeEventsCountMap[uNode];

		unordered_map<ui, ui> nodeEdgeCounts =  nodeNodeCountIterator->second;

		unordered_map<ui, ui>::iterator nodeEdgeCountsIterator;

		for(nodeEdgeCountsIterator = nodeEdgeCounts.begin(); nodeEdgeCountsIterator != nodeEdgeCounts.end(); nodeEdgeCountsIterator++)
		{
			ui vNode = nodeEdgeCountsIterator->first;
			ui Nuv = nodeEdgeCountsIterator->second;

			double wuv = ((Nuv + 0.01) * 1.0)/(Nu + 1);
			userUserInfluence[uNode][vNode] = wuv;

			avgVal += wuv;
			count += 1;
		}
		// nodeNodeCount[uNode].clear();
	}

	nodeNodeCount.clear();
	cout << "Avg Wuv value = " << (avgVal * 1.0) / count  << " Over Edges (edges with transaction) = " << count << endl;

	return 0;
}




int updateNodeNodeCountMap()
{
	li eventNode, parentEvent, parentNode;
	
	nodeNodeCount.clear();

	// depending on the current parent assignment we get the nodeNode Count... i.e. Nuv values... for each edge...
	for(unsigned int i = 0; i < newSyntheticEvents.size(); i++)
	{
		eventNode = newSyntheticEvents[i][1];
		parentEvent = newSyntheticEvents[i][2];
		if(parentEvent >= 0)
		{
			parentNode = newSyntheticEvents[parentEvent][1];
			nodeNodeCount[parentNode][eventNode]++;
		}
	}

	return 0;
}


///////////////////////////////
///////  SAMPLE TOPIC    //////
///////////////////////////////

int sampleTopicAssignment(int ITE)
{
	cout << "Sample Topic Assignment\n";

	for(unsigned int i = 0; i < newSyntheticEvents.size(); i++)	
	{
		li eventIndex = newSyntheticEvents[i][0];
		// li eventIndex = i;
		li eventNode = newSyntheticEvents[i][1];
		li eventParent = newSyntheticEvents[i][2];			// should be -2 for the first iteration...
		li eventTopic = newSyntheticEvents[i][3];				// should be -2 for the first iteration...	

		vector<ui> doc = allDocs[i];

		li assignedTopic;
		
		decreamentCountFromMatrices(eventIndex, eventNode, eventParent, eventTopic, doc, 1);
		assignedTopic =  getSampledTopicAssignment(eventIndex, eventNode, eventParent, doc, ITE);
		increamentCountToMatrices(eventIndex, eventNode, eventParent, assignedTopic, doc, 1);

		newSyntheticEvents[eventIndex][3] = assignedTopic;

		// getting every 10th iteration topic assigned, i.e. the topic assigned at each Markov Chain... 
		if(ITE > BURN_IN && ITE%10 == 0)
		{
			topicDistEvery10thIter[eventIndex].push_back(assignedTopic);
		}
	}
	// to get some stats...
	printTopicTopicCount();

	return 0;
}

int printTopicTopicCount()
{
	int sumAll = 0;
	for(ui i = 0; i < numTopics; i++)
	{
		for(ui j = 0; j < numTopics; j++)
		{
			// cout << NTTCountVec[i][j] << " ";
			sumAll += NTTCountVec[i][j];
		}
		// cout << "\n";
	}
	cout << "Sum Over all i-j -- " << sumAll << "\n";
	
	if(sumAll > countInteractions)
	{
		cout << "Sum all is greater than Count Interactions (For first iteration it can be...)\n";
		cout << sumAll << " " << countInteractions << "\n";
		// exit(0);
	}

	return 0;
}


li getSampledTopicAssignment(li eventIndex, li eventNode, li eventParent, vector<ui> doc, int ITE)
{
	li assignedTopic = -1;

	vector <double> calculatedProbVec (numTopics, 0.0);
	
	// we need the hist of topics assigned to the child events... For the second/middle term specifically...
	unordered_map <int, ui> childEventTopicsHist = getHistOfTopicsOverChildEvents(eventIndex);

	// calculating probability for each possible topic
	for(li topic = 0; topic < numTopics; topic++)
	{
		// cout << "Topic iterator at - " << topic << endl;
		double firstTerm = 0.0;
		double middleTerm = 0.0;
		double thirdTerm = 0.0;

		// cout << "eventIndex -- " << eventIndex << endl;
		// first term -- related to the parents topic... or else to the node Topic preference in case of no parent...
		firstTerm = getFirstTermOfTopicAssignmentCondProb(eventNode, eventParent, topic);
		
		if(childEventTopicsHist.size() > 0)
		{
			middleTerm = getMiddleTermOfTopicAssignmentCondProb(topic, childEventTopicsHist, eventIndex, eventParent);
		}
		else
		{
			// if there is no child event... then there is no event to iterate on... 
			// and the probability should be calculated just on the basis of the parent event and words in the event...
			middleTerm = 0;
		}
		// middleTerm = 0;
		// hist of words for each doc do not change over iteration...  
		vector<ui> wordHistVec = wordHistAllDocsVector[eventIndex];
		thirdTerm = getThirdTermOfTopicAssignmentCondProb(topic, doc, wordHistVec);
		// thirdTerm = 0;

		calculatedProbVec[topic] = firstTerm + middleTerm + log(thirdTerm);
		// cout << "Topic iterator at done -- " << topic << endl;
	}

	assignedTopic = getSampleFromMultinomial(calculatedProbVec);

	if(ITE > BURN_IN && ITE%10 == 0)
	{
		for(ui i = 0; i < avgTopicProbVector[eventIndex].size(); i++)
		{
			avgTopicProbVector[eventIndex][i] += calculatedProbVec[i];
		}
	}

	calculatedProbVec.clear();
	// cout << assignedTopic << "\n";

	return assignedTopic;
}

double getFirstTermOfTopicAssignmentCondProb(li eventNode, li eventParent, li topic)
{
	// not calculating the denominator as normalization is going to be same for the all the topics...
	// irrespective of whether the event has a parent or not...
	double finalFirstTerm = 0;

	if(eventParent >= 0)
	{
		li eventParentTopic = newSyntheticEvents[eventParent][3];
		// finalFirstTerm = hyperBeta[topic] + NTTCountVec[eventParentTopic][topic];
		finalFirstTerm = hyperBeta[topic] + NTTCountVec[eventParentTopic][topic];
	}
	else if(eventParent == -1)
	{
		// no parent
		finalFirstTerm = hyperGamma[topic] + NUTCountVec[eventNode][topic];
		// cout << eventNode << " " << topic << endl;
	}
	else
	{
		// This would not get executed if there is some parent initialization...
		// if the parent is invalid or not assigned yet...
		// finalFirstTerm = log(hyperGamma[topic]);
		finalFirstTerm = hyperGamma[topic];
	}

	double logVal = log(finalFirstTerm);
	// cout << logVal << endl;
	return logVal;
}


double getMiddleTermOfTopicAssignmentCondProb(li topic, unordered_map <int, ui> childEventTopicsHist, li eventIndex, li eventParent)
{
	// not calculating this middle term distinctly when parentTopics = currentTopic
	// as there is no much change... the equations in both the case are almost similar... 
	double finalMiddleValue = 0.0;
	double finalMiddleValueNume = 0.0;
	double finalMiddleValueDenom = 0.0;

	li eventParentTopic;

	if(eventParent >= 0)
	{
		eventParentTopic = newSyntheticEvents[eventParent][3];
	}
	else
	{
		eventParentTopic = -2;				// some invalid value...
	}

	unordered_map <int, ui>::iterator childEventTopicsHistIt;

	// int childEventTopicsHistCount = 0;

	for(childEventTopicsHistIt = childEventTopicsHist.begin(); childEventTopicsHistIt != childEventTopicsHist.end(); childEventTopicsHistIt++)
	{
		// topic lPrime in the child event...
		li lPrime = childEventTopicsHistIt->first;

		// if there is some topic assigned to the child events...  
		if(lPrime > -1)
		{
			// count of the topics in the child events...
			li lPrimeCountInChilds = childEventTopicsHistIt->second;

			// For a particular possible topic the base term remains same 
			// over the count of the topics in the neighborhood/child events...

			// double baseTerm = betaLPrime + ttCountWithoutChildEvents;
			double baseTerm = hyperBeta[lPrime] + NTTCountVec[topic][lPrime];

			// value for each topic...
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
	double baseTermDenom = sumBeta + NTTSumTopicsVec[topic];

	ui numChildEvents = 0;

	vector <ui> cEvents;

	try
	{
		// get the number of child events from the updated child events map...
		numChildEvents = childEventsMap[eventIndex].size();
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
unordered_map <int, ui> getHistOfTopicsOverChildEvents(li eventIndex)
{
	unordered_map <int, ui> topicCounts;

	vector <ui> childEventsList;
	try
	{
		childEventsList =  childEventsMap.at(eventIndex);
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
			li childEventTopic = newSyntheticEvents[childEventIndex][3];

			topicCounts[childEventTopic] += 1;
		}
	}

	return topicCounts;
}  


double getThirdTermOfTopicAssignmentCondProb(li topic, vector<ui> doc, vector<ui> wordHistVec)
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
		double baseterm = alphaWord + NTWCountVec[topic][word];

		long double valueForWord = 1.0;

		for(ui j = 0; j < wordCount; j++)
		{
			valueForWord *= (baseterm + j);
		}

		finalThirdValueNume *= valueForWord;
		i++;
	}

	double baseTermDenom = sumAlpha + NTWSumWordsVec[topic];

	for(ui i = 0; i < doc.size(); i++)
	{
		finalThirdValueDenom *= (baseTermDenom + i);
	}

	finalThirdValue = finalThirdValueNume / finalThirdValueDenom;

	return finalThirdValue;
}

// this is required for calculating the lastterm of the probability for topic assignment...
int createWordHistForAllDocs()
{
	for(unsigned int i = 0; i < allDocs.size(); i++)
	{
		vector<ui> doc = allDocs[i];

		unordered_map<ui, ui> wordHist;
		unordered_map<ui, ui>::iterator wordHistIt;
		wordHist = getHistOfWordsOverWordsFromDoc(doc);

		// wordHistAllDocs[i] = wordHist;
		vector <ui> tempWordCounts;
		for(wordHistIt = wordHist.begin(); wordHistIt != wordHist.end(); wordHistIt++)
		{
			tempWordCounts.push_back(wordHistIt->first);
			tempWordCounts.push_back(wordHistIt->second);
		}
		wordHistAllDocsVector.push_back(tempWordCounts);
	}

	return 0;
}

// this is required for calculating the lastterm of the probability for topic assignment...
unordered_map <ui, ui> getHistOfWordsOverWordsFromDoc(vector<ui> doc)
{
	unordered_map <ui, ui> wordCounts;

	for(unsigned int i = 0; i < doc.size(); i++)
	{
		li word = doc[i];

		wordCounts[word] += 1;
	}

	return wordCounts;
}

///////////////////////////////
///////  SAMPLE PARENT   //////
///////////////////////////////

int sampleParentAssignment(int ITE)
{
	cout << "Sampling Parent Assignment ... \n";

	// not using nodeNodeCombUpdateInfluence for now...
	// this is to keep the edges that have change in the parent...
	nodeNodeCombUpdateInfluence.clear();
	cout << "nodeNodeUpdates cleared count = " << nodeNodeCombUpdateInfluence.size() << "\n";
	// logLikelihood = 0;

	li previousParentNode, assignedParentNode;

	for(unsigned int i = 0; i < newSyntheticEvents.size(); i++)
	{
		li eventIndex = newSyntheticEvents[i][0];
		li eventNode = newSyntheticEvents[i][1];

		li eventParent = newSyntheticEvents[i][2];
		li eventTopic = newSyntheticEvents[i][3];

		double eventTime = eventIndexTimestamps[i];

		vector<ui> doc;
		// cout << "eventIndex -- " << eventIndex << endl;
		// this is to understand/update the change in the parent assignment...
		li previousParent = eventParent;

		if(previousParent >= 0)
		{
			previousParentNode = newSyntheticEvents[previousParent][1];
		}
		else
		{
			previousParentNode = eventNode;
		}

		decreamentCountFromMatrices(eventIndex, eventNode, eventParent, eventTopic, doc, 0);
		li assignedParent = getSampledParentAssignment(eventTime, eventNode, eventIndex, eventTopic, ITE);
		increamentCountToMatrices(eventIndex, eventNode, assignedParent, eventTopic, doc, 0);
		// assign the new sampled parent event to the event...
		
		newSyntheticEvents[eventIndex][2] = assignedParent;

		if(assignedParent >= 0)
		{
			assignedParentNode = newSyntheticEvents[assignedParent][1];
		}
		else
		{
			assignedParentNode = eventNode;
		}

		if(previousParent == -1 && assignedParent >= 0)
		{
			nodeNodeCombUpdateInfluence[assignedParentNode].push_back(eventNode);
		}
		else if(previousParent >= 0 && assignedParent == -1)
		{
			nodeNodeCombUpdateInfluence[previousParentNode].push_back(eventNode);
		}
		else if(previousParent >= 0 && assignedParent >= 0)
		{
			nodeNodeCombUpdateInfluence[assignedParentNode].push_back(eventNode);
			nodeNodeCombUpdateInfluence[previousParentNode].push_back(eventNode);
		}

		// if(i % 10000 == 0)
			// printf("Parent Assignment -- %d\n", i);

		if(ITE > BURN_IN && ITE%10 == 0)
		{
			parentDistEvery10thIter[eventIndex].push_back(assignedParent);
		}
	}

	// cout << "\nLog Likelihood = " << logLikelihood << "\n";

	return 0;
}

li getSampledParentAssignment(double eventTime, li eventNode, li eventIndex, li eventTopic, int ITE)
{
	li assignedParent = -2;

	vector <ui> possibleParentEvents = allPossibleParentEvents[eventIndex];
	vector<double> possibleParentExp = allPossibleParentEventsExponentials[eventIndex];
	
	vector <double> calculatedProbVec;

	if(possibleParentEvents.size() > 0)
	{
		calculatedProbVec = populateCalculatedProbVec(possibleParentEvents, possibleParentExp, eventNode, eventTopic, eventTime, ITE);

		// cout << "EventIndex -- " << eventIndex << "-- calculatedProbVec.size -- " << calculatedProbVec.size() << endl;

		// ui sampledIndex = getSampleFromMultinomial(calculatedProbVec);
		ui sampledIndex = getSampleFromDiscreteDist(calculatedProbVec);

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
		for(ui i = 0; i < avgProbParForAllEvents[eventIndex].size(); i++)
		{
			avgProbParForAllEvents[eventIndex][i] += calculatedProbVec[i];
		}
	}

	if(assignedParent > (int)allEvents.size()-1)
	{
		cout << "some issue with parent sampling...\n";
		exit(0);
	}

	return assignedParent;
}

vector<double> populateCalculatedProbVec(vector <ui> possibleParentEvents, vector<double> possibleParentExp, li eventNode, li eventTopic, double eventTime, int ITE)
{

	// the last index is for probability of self event being the parent...
	vector <double> calculatedProbVec(possibleParentEvents.size() + 1, 0.0);
	
	// double firstTermNume, firstTermDenom, firstTerm, secondTerm;
	double firstTerm, secondTerm;

	double normalizationFactor = 0;

	unsigned int k;
	for(k = 0; k < possibleParentEvents.size(); k++)
	{
		vector <li> possParentEvent = newSyntheticEvents[possibleParentEvents[k]];
		// li possParentIndex = possParentEvent[0];
		li possParentNode = possParentEvent[1];
		li possParentEventTopic = possParentEvent[3];

		// this is the parent dependent term...
		firstTerm = getFirstTermOfParentAssignment(possParentEventTopic, eventTopic);

		double Wuv;

		Wuv = userUserInfluence[possParentNode][eventNode];

		if(Wuv > 0)
		{
			double temporalDecay = possibleParentExp[k];
			// this is Wuv times the temporal decay term...
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
	firstTerm = getFirstTermOfParentAssignmentNoParent(eventNode, eventTopic);
	
	if(userBaseRateMap.find(eventNode) != userBaseRateMap.end())
	{
		secondTerm = userBaseRateMap[eventNode];
	}
	else
	{
		// secondTerm = defaultMuVal;
		secondTerm = 0;
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
	// get the normalized probability vector...
	for(ui i = 0; i < calculatedProbVec.size(); i++)
	{
		calculatedProbVec[i] = calculatedProbVec[i] / normalizationFactor;
	}

	return calculatedProbVec;
}

double getFirstTermOfParentAssignment(li possParentEventTopic, li eventTopic)
{
	double firstTerm = 0.0;
	double firstTermNume = 0.0;
	double firstTermDenom = 1.0;
	
	firstTermNume = hyperBeta[eventTopic] + NTTCountVec[possParentEventTopic][eventTopic];
	
	firstTermDenom = sumBeta + NTTSumTopicsVec[possParentEventTopic];

	firstTerm = firstTermNume / firstTermDenom;

	return firstTerm;
}

double getFirstTermOfParentAssignmentNoParent(li eventNode, li eventTopic)
{
	double firstTerm = 0.0;
	double firstTermNume = 0.0;
	double firstTermDenom = 1.0;

	firstTermNume = hyperGamma[eventTopic] + NUTCountVec[eventNode][eventTopic];
	firstTermDenom = sumGamma + NUTSumTopicsVec[eventNode];
	firstTerm = firstTermNume / firstTermDenom;

	return firstTerm;
}

vector<ui> getPossibleParentEvents(li eventNode, li eventIndex)
{
	vector <ui> possibleParentEvents = allPossibleParentEvents[eventIndex];  
	return possibleParentEvents;
}


// Decreament Counters
int decreamentCountFromMatrices(li eventIndex, li eventNode, li eventParent, li eventTopic, vector <ui> doc, bool topicSampling)
{
	// cout << "decreament count from topic-topic matrix.. from the cell parentTopic -> eventTopic";
	// decreament only if the eventTopic is valid...
	if(eventTopic >= 0)
	{
		int eventParentTopic;

		if(eventParent >= 0)
		{
			eventParentTopic = newSyntheticEvents[eventParent][3];
			// printUnorderedMap("ttopic");
			decreamentCountFromTopicTopic(eventParentTopic, eventTopic);
			// printUnorderedMap("ttopic");
		}
		else if(eventParent == -1)
		{
			// printUnorderedMap("utopic");
			decreamentCountFromUserTopic(eventNode, eventTopic);
			// printUnorderedMap("utopic");
		}

		if(topicSampling == 1)
		{
			// decreament the count of words corresponding to this event...
			decreamentCountFromTopicWord(eventTopic, doc);
			
			// decreament counts corresponding to child events for topic sampling...
			decreamentCountsFromChildEvents(eventTopic, eventIndex);
		}
	}

	return 0;
}


int decreamentCountFromTopicTopic(li eventParentTopic, li eventTopic)
{
	// cout << "decreament count from topic-topic matrix.. from the cell parentTopic -> eventTopic";

	if(NTTCountVec[eventParentTopic][eventTopic] > 0)
	{
		NTTCountVec[eventParentTopic][eventTopic] -= 1;
	}

	if(NTTSumTopicsVec[eventParentTopic] > 0)
	{
		NTTSumTopicsVec[eventParentTopic] -= 1;
	}

	return 0;
}


int decreamentCountFromUserTopic(li eventNode, li eventTopic)
{
	// cout << "decreament count from user-topic matrix.. from the cell user -> eventTopic";

	if(NUTCountVec[eventNode][eventTopic] > 0)
	{
		NUTCountVec[eventNode][eventTopic] -= 1;	
	}

	if(NUTSumTopicsVec[eventNode] > 0)
	{
		NUTSumTopicsVec[eventNode] -= 1;	
	}

	return 0;
}

int decreamentCountFromTopicWord(li eventTopic, vector <ui> doc)
{
	if(NTWSumWordsVec[eventTopic] >= doc.size())
	{
		NTWSumWordsVec[eventTopic] -= doc.size();
	}
	else
	{
		// cout << "Somethings wrong with the Topic Word Counts... \n";
		NTWSumWordsVec[eventTopic] = 0;
	}

	for(unsigned int i = 0; i < doc.size(); i++)
	{
		if(NTWCountVec[eventTopic][doc[i]] > 0)
		{
			NTWCountVec[eventTopic][doc[i]] -= 1;
		}
	}

	return 0;
}


int decreamentCountsFromChildEvents(li eventTopic, li eventIndex)
{
	// vector <li> childEventsList;
	vector <ui> childEventsList;
	try
	{
		childEventsList =  childEventsMap.at(eventIndex);
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

			li childEventTopic = newSyntheticEvents[childEventIndex][3];

			if(childEventTopic != -1)
			{
				decreamentCountFromTopicTopic(eventTopic, childEventTopic);
			}
		}
	}

	return 0;
}


// Increament Counters
int increamentCountToMatrices(li eventIndex, li eventNode, li eventParent, li eventTopic, vector<ui> doc, bool topicSampling)
{
	// cout << "increament count to topic-topic matrix... in the cell parentTopic -> assignedTopic";

	if(eventParent >= 0)
	{
		int eventParentTopic = newSyntheticEvents[eventParent][3];
		// printUnorderedMap("ttopic");
		// cout << "increamenting...\n";
		increamentCountInTopicTopic(eventParentTopic, eventTopic);
		// printUnorderedMap("ttopic");
	}
	else if (eventParent == -1)
	{
		// printUnorderedMap("utopic");
		// cout << "increamenting...\n";
		increamentCountInUserTopic(eventNode, eventTopic);
		// printUnorderedMap("utopic");
	}

	if(topicSampling == 1)
	{
		// increament counts in topic-word...
		increamentCountInTopicWord(eventTopic, doc);
		// printUnorderedMap("psitopic");
		
		// increament counts corresponding to child events...
		increamentCountsForChildEvents(eventTopic, eventIndex);
	}

	return 0;
}

int increamentCountInTopicTopic(li eventParentTopic, li eventTopic)
{

	NTTSumTopicsVec[eventParentTopic] += 1;
	NTTCountVec[eventParentTopic][eventTopic] += 1;

	if(NTTCountVec[eventParentTopic][eventTopic] > newSyntheticEvents.size())
	{
		cout << "Somethings wrong.. the size of NTT is more than number of events...\n";
		cout << NTTCountVec[eventParentTopic][eventTopic] << " " << newSyntheticEvents.size() << "\n";
		exit(0);
	}


	if(NTTSumTopicsVec[eventParentTopic] > newSyntheticEvents.size())
	{
		cout << "Somethings wrong.. the size of NTT (sum) is more than number of events...\n";
		cout << NTTSumTopicsVec[eventParentTopic] << " " << newSyntheticEvents.size() << "\n";
		exit(0);
	}

	// NTSumTopics[eventParentTopic] += 1;

	// li cellIndex = getCellKey(eventParentTopic, eventTopic);
	// NTTCountMatrix[cellIndex] += 1;

	return 0;
}

int increamentCountInUserTopic(li eventNode, li eventTopic)
{
	NUTSumTopicsVec[eventNode] += 1;
	NUTCountVec[eventNode][eventTopic] += 1;
	
	if(NUTCountVec[eventNode][eventTopic] > newSyntheticEvents.size())
	{
		cout << "Somethings wrong.. the size of NUT is more than number of events...\n";
		cout << NUTCountVec[eventNode][eventTopic] << " " << newSyntheticEvents.size() << "\n";
		exit(0);
	}


	if(NUTSumTopicsVec[eventNode] > newSyntheticEvents.size())
	{
		cout << "Somethings wrong.. the size of NUT (sum) is more than number of events...\n";
		cout << NUTSumTopicsVec[eventNode] << " " << newSyntheticEvents.size() << "\n";
		exit(0);
	}


	// NUSumTopics[eventNode] += 1;

	// li cellIndex = getCellKey(eventNode, eventTopic);
	// NUTCountMatrix[cellIndex] += 1;

	return 0;
}

int increamentCountInTopicWord(li eventTopic, vector<ui> doc)
{
	NTWSumWordsVec[eventTopic] += doc.size();
	// NTSumWords[eventTopic] += doc.size();

	if(NTWSumWordsVec[eventTopic] > totalWords)
	{
		cout << "Total words = " << totalWords << " NTW sum = " << NTWSumWordsVec[eventTopic] << "\n";
		exit(0);
	}

	for(unsigned int i = 0; i < doc.size(); i++)
	{
		NTWCountVec[eventTopic][doc[i]] += 1;

		if(NTWCountVec[eventTopic][doc[i]] > totalWords)
		{
			cout << "Total words = " << totalWords << " NTW = " << NTWCountVec[eventTopic][doc[i]] << " " << doc[i] << "\n";
			exit(0);
		}
	}
	
	return 0;
}

int increamentCountsForChildEvents(li eventTopic, li eventIndex)
{
	// vector <li> childEventsList;
	vector <ui> childEventsList;
	try
	{
		childEventsList = childEventsMap.at(eventIndex);
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

			li childEventTopic = newSyntheticEvents[childEventIndex][3];

			if(childEventTopic != -1)
			{
				increamentCountInTopicTopic(eventTopic, childEventTopic);
			}
		}
	}

	return 0;
}

///////////////////////////////
//////// UTIL FUNCTIONS ///////
///////////////////////////////

li getCellKey(li firstNum, li secondNum)
{
	// string keyString;

	// string firstUnderscore = to_string(firstPart) + "_";
	// keyString = firstUnderscore + to_string(secondPart);

	// using Szudziks Function....
	// a >= b ? a * a + a + b : a + b * b
	li keyVal = firstNum >= secondNum ? firstNum * firstNum + firstNum + secondNum : firstNum + secondNum * secondNum;



	// return keyString;
	return keyVal;
}

int getSampleFromMultinomial(vector<double> calculatedProbVec)
{
	int assignedInd;

	// get normalized prob vector...
	vector <double> normalizedProbVector = getNormalizedLogProb(calculatedProbVec);
	
	assignedInd = getSampleFromDiscreteDist(normalizedProbVector);

	// get the log probability of the assigned topic and add it to the loglikelihood... 
	logLikelihood += log(normalizedProbVector[assignedInd]);

    return assignedInd;
}

/*
// This is prone to overflow....
vector<double> getNormalizedLogProb(vector<double> calculatedProbVec)
{
	vector<double> normalizedProbVector(calculatedProbVec.size(), 0.0);

	// each of the terms will be normalized as:
	// log(pi) - ( log(mi) + log(sum (exp (log(pi) - log(mi)))) )

	// double maxTerm = *max_element(calculatedProbVec.begin(), calculatedProbVec.end());

	double sumOverExp = 0.0;
	// sum(exp(log(pi) - log(mi)))
	for(unsigned int i = 0; i < calculatedProbVec.size(); i++)
	{
		// sumOverExp += exp(calculatedProbVec[i] - maxTerm);
		sumOverExp += exp(calculatedProbVec[i]);
	}

	for(unsigned int i = 0; i < calculatedProbVec.size(); i++)
	{
		// normalizedProbVector[i] = calculatedProbVec[i] - (log(maxTerm) + log(sumOverExp));
		// normalizedProbVector[i] = calculatedProbVec[i] - (log(sumOverExp));
		normalizedProbVector[i] = exp(calculatedProbVec[i] - (log(sumOverExp)));
	}

	return normalizedProbVector;
}
*/

vector<double> getNormalizedLogProb(vector<double> calculatedProbVec)
{
	// use the log sum exp trick here to avoid underflow...
	vector<double> normalizedProbVector(calculatedProbVec.size(), 0.0);

	// each of the terms will be normalized as:
	// log(pi) - ( log(mi) + log(sum (exp (log(pi) - log(mi)))) )

	double maxTerm = *max_element(calculatedProbVec.begin(), calculatedProbVec.end());

	double sumOverExp = 0.0;
	// sum(exp(log(pi) - log(mi)))
	for(unsigned int i = 0; i < calculatedProbVec.size(); i++)
	{
		sumOverExp += exp(calculatedProbVec[i] - maxTerm);
	}

	for(unsigned int i = 0; i < calculatedProbVec.size(); i++)
	{
		normalizedProbVector[i] = exp(calculatedProbVec[i] - (maxTerm + log(sumOverExp)));
		// normalizedProbVector[i] = exp(calculatedProbVec[i] - (log(sumOverExp)));
	}

	return normalizedProbVector;
}

int getSampleFromDiscreteDist(vector<double> normalizedProbVector)
{
	// remove .. later
	random_device rd;
	mt19937 gen(rd());
	
	// default_random_engine gen;

	discrete_distribution<int> distribution(normalizedProbVector.begin(), normalizedProbVector.end());
	int ind = distribution(gen);

	return ind;
}


double getSampleFromGamma(double alpha, double beta)
{
	random_device rd;
	mt19937 gen(rd());
	gamma_distribution<double> distribution(alpha, beta);

	double gammaVal = distribution(gen);

	return gammaVal;
}


///////////////////////////////
////////   VALIDATIONS   ///////
///////////////////////////////


int writeEvery10ItersTopicAssignmentToFile()
{
	ofstream topicAssignmentsFile;
	topicAssignmentsFile.open("B5_topicAssignments_OurModel_estAll.txt");

	cout << "Writing topic assingments to file...\n";

	map<int, vector <int> >::iterator topicDistEvery10thIterIterator;

	for(topicDistEvery10thIterIterator = topicDistEvery10thIter.begin(); topicDistEvery10thIterIterator != topicDistEvery10thIter.end(); topicDistEvery10thIterIterator++)
	{
		int eventId = topicDistEvery10thIterIterator->first;
		vector<int> assignments = topicDistEvery10thIterIterator->second;

		topicAssignmentsFile << eventId << " ";

		for(ui i = 0; i < assignments.size(); i++)
		{
			topicAssignmentsFile << assignments[i] << " ";			
		}

		topicAssignmentsFile << endl;
	}

	return 0;
}

int writeEvery10ItersParentAssignmentToFile()
{
	ofstream parentAssignmentsFile;
	parentAssignmentsFile.open("B5_parentAssignments_OurModel_estAll.txt");

	cout << "Writing parent assignments to file ...\n";

	map<int, vector <int> >::iterator parentDistEvery10thIterIterator;

	for(parentDistEvery10thIterIterator = parentDistEvery10thIter.begin(); parentDistEvery10thIterIterator != parentDistEvery10thIter.end(); parentDistEvery10thIterIterator++)
	{
		int eventId = parentDistEvery10thIterIterator->first;
		vector<int> assignments = parentDistEvery10thIterIterator->second;

		parentAssignmentsFile << eventId << " ";

		for(ui i = 0; i < assignments.size(); i++)
		{
			parentAssignmentsFile << assignments[i] << " ";			
		}

		parentAssignmentsFile << endl;
	}

	return 0;
}

int writeAvgProbVectorsToFile()
{
	ofstream avgProbVecFile;
	avgProbVecFile.open("B5_avgParProb_OurModel_estAll.txt");

	cout << "Writing Avg Prob Vectors to file ...\n";
	// format
	// eid numofProbValues pid1 probval1 pid2 probval2 ...
	int takeAvgOver = (ITERATIONS - BURN_IN - 1)/10;
	for(ui i = 0; i < avgProbParForAllEvents.size(); i++)
	{
		avgProbVecFile << avgProbParForAllEvents[i].size() << " " << i << " ";

		ui j = 0;
		for(j = 0; j < avgProbParForAllEvents[i].size() - 1; j++)
		{
			avgProbVecFile << allPossibleParentEvents[i][j] << " " <<  avgProbParForAllEvents[i][j] / takeAvgOver << " ";
		}

		avgProbVecFile << i << " " << avgProbParForAllEvents[i][j] / takeAvgOver << "\n";
	}

	avgProbVecFile.close();
	return 0;
}


int writeAvgTopicProbVectorsToFile()
{
	ofstream avgTopicProbVecFile;
	// avgProbVecFile.open("avgParProb_OurModel.txt");
	avgTopicProbVecFile.open("B5_avgTopicProb_OurModel_estAll.txt");

	cout << "Writing Avg Prob Vectors to file ...\n";
	// format
	// eid numofProbValues pid1 probval1 pid2 probval2 ...
	for(ui i = 0; i < avgTopicProbVector.size(); i++)
	{
		avgTopicProbVecFile << avgTopicProbVector[i].size() << " " << i << " ";

		ui j = 0;
		for(j = 0; j < avgTopicProbVector[i].size(); j++)
		{
			ui topicId = j;
			avgTopicProbVecFile << topicId << " " <<  avgTopicProbVector[i][j] / 10 << " ";
		}

		avgTopicProbVecFile << "\n";
	}

	avgTopicProbVecFile.close();
	return 0;
}


int writeTopicTopicInteractionToFile()
{
	ofstream topicTopicIntOutFile;
	topicTopicIntOutFile.open("B5_topicTopicInteraction_OurModel_estAll.txt");

	cout << "Writing topic-topic interactions ...\n";

	for(ui i = 0; i < NTTCountVec.size(); i++)
	{
		topicTopicIntOutFile << i << " ";

		for(ui j = 0; j < NTTCountVec[i].size(); j++)
		{
			topicTopicIntOutFile << NTTCountVec[i][j] << " ";
		}

		topicTopicIntOutFile << "\n";
	}

	topicTopicIntOutFile.close();
	
	return 0;
}

int writeUserUserInfluenceToFile()
{
	ofstream uuinfFile;
	uuinfFile.open("B5_userUserInf_OurModel_estAll.txt");

	cout << "Writing user user influences...\n";

	unordered_map <ui, unordered_map<ui, double> >::iterator userUserInfluenceIterator;

	for(userUserInfluenceIterator = userUserInfluence.begin(); userUserInfluenceIterator != userUserInfluence.end(); userUserInfluenceIterator++)
	{
		ui uNode = userUserInfluenceIterator->first;
		unordered_map<ui, double> userInf = userUserInfluenceIterator->second;

		if(userInf.size() > 0)
		{
			uuinfFile << userInf.size() << " " << uNode << " ";

			unordered_map<ui, double>::iterator userInfIterator;

			for(userInfIterator = userInf.begin(); userInfIterator != userInf.end(); userInfIterator++)
			{
				uuinfFile << userInfIterator->first << " " << userInfIterator->second << " "; 
			}
			
			uuinfFile << "\n";
		}

	}

	uuinfFile.close();

	return 0;
}


int writeUserBaseRatesToFile()
{

	ofstream ubrFile;
	ubrFile.open("B5_userBaseRate_OurModel.txt");

	map<ui, double>::iterator userBaseRateMapIt;

	for(userBaseRateMapIt = userBaseRateMap.begin(); userBaseRateMapIt != userBaseRateMap.end(); userBaseRateMapIt++)
	{
		int uid = userBaseRateMapIt->first;
		double ubrUid = userBaseRateMapIt->second;

		// ubrFile << uid << " " << ubrUid << "\n";
		ubrFile << setprecision(10) << uid << " " << ubrUid << "\n";
	}

	ubrFile.close();

	return 0;
}


int writeNodeNodeCountAndNodeCount()
{
	ofstream wuvAnalyze;
	wuvAnalyze.open("B5_Nuv_Nu_Wuv_OurModel.txt");

	unordered_map<ui, unordered_map<ui, ui> >::iterator nodeNodeCountIt;

	for (nodeNodeCountIt = nodeNodeCount.begin(); nodeNodeCountIt != nodeNodeCount.end(); nodeNodeCountIt++)
	{
		ui uNode = nodeNodeCountIt->first;

		unordered_map<ui, ui> edgeCountOverNode = nodeNodeCountIt->second;

		unordered_map<ui, ui>::iterator edgeCountOverNodeIt;

		for(edgeCountOverNodeIt = edgeCountOverNode.begin(); edgeCountOverNodeIt != edgeCountOverNode.end(); edgeCountOverNodeIt++)
		{
			ui vNode = edgeCountOverNodeIt->first;
			ui Nuv = edgeCountOverNodeIt->second;

			ui Nu = nodeEventsCountMap[uNode];

			double wuv = ((Nuv + 0.01) * 1.0)/(Nu + 1);


			wuvAnalyze << vNode << " " << uNode << " " << Nuv << " " << Nu << " " << wuv << "\n";
		}
	}

	wuvAnalyze.close();

	return 0;
}

///////////////////////////////
///////// READING DATA ////////
//////// AND PARAMETERS ///////
//////// FROM FILES ///////////
///////////////////////////////

// vector < vector <double> > getSyntheticEventsFromFile(string fileName)
vector < vector <li> > getSyntheticEventsFromFile(string fileName)
{
	vector < vector <li> > localAllEvents;

	ifstream allEventsFile;
    allEventsFile.open(fileName);

    string line;
    stringstream ss;

    unsigned i = 0;

    maxLevel = 0;

    if(allEventsFile.is_open())
    {
    	while(getline(allEventsFile, line))
		{
			double eventTime;
			li eventNode, eventParent, eventTopic;
			ui level;

			ss.clear();
			ss.str("");

			ss << line;
			ss >> eventTime >> eventNode >> eventParent >> eventTopic >> level;

			eventIndexTimestamps[i] = eventTime;

			vector<li> tempEvent;
			tempEvent.push_back(i);						// eventIndex...
			tempEvent.push_back(eventNode);
			tempEvent.push_back(eventParent);			
			tempEvent.push_back(eventTopic);

			levelInfo.push_back(level);
			
			if (maxLevel < level)
			{
				maxLevel = level;
			}

			localAllEvents.push_back(tempEvent);

			line.clear();

			i++;

			tempEvent.clear();
		}
    }
	else
	{
		cout << "Error opening file -- " << fileName << "\n";
	}
	cout << "read " << localAllEvents.size() << " events\n";

	cout << "First event time -- " << eventIndexTimestamps[0] << "\n";

	allEventsFile.close();
	return localAllEvents;
}

vector < vector <ui> > getSyntheticDocsFromFile(string fileName)
{
	vector < vector <ui> > allDocs;

	double avgDocLength = 0;

	ifstream allDocsFile;
    allDocsFile.open(fileName);

    string line;
    stringstream ss;
    li i = 0;

    if(allDocsFile.is_open())
    {
	    while(getline(allDocsFile, line))
		{
			ss.clear();
			ss.str("");

			ss << line;

			int word;
			vector <ui> currDoc;
			while(ss >> word)
			{
				// vector<double> tempEvent;
				currDoc.push_back(word);
			}
			
			eventIndexDocSizeMap[i] = currDoc.size();
			i++;

			allDocs.push_back(currDoc);
			line.clear();

			avgDocLength += currDoc.size();

			// if(i > 100000)
				// break;
		}
		allDocsFile.close();

		totalWords = avgDocLength;
		cout << "Total Number of words = " << totalWords << "\n";

		avgDocLength = avgDocLength/allDocs.size();
		
		cout << "Read " << allDocs.size() << " documents\n";
		cout << "Avg Doc Length = " << avgDocLength << endl;
    }
    else
    {
    	cout << "Error opening file -- " << fileName << "\n";
    }
    cout << "read " << allDocs.size() << " docs\n";
	return allDocs;
}

// unordered_map <ll, vector <ll> > readIntVectorMapFromFile(string fileName)
vector < vector <ui> > readIntVectorMapFromFile(string fileName)
{
	ifstream intVectorFile;
	intVectorFile.open(fileName);
	// intVectorFile.open("mapped_users_followers_graph_sorted_1000000_users.txt");

	// unordered_map <ll, vector <ll> > intVectorMap;
	vector < vector <ui> > intVectorMap(maxNumNodes);

	string line;
	stringstream ss;
	ui nodeId;
	ui follower;
	ui count; 
	ui mapsize = 0;

	if (intVectorFile.is_open())
	{
		while(getline(intVectorFile, line))
		{
			ss.clear();
			ss.str("");

			ss << line;
			ss >> count >> nodeId;

			vector <ui> tempVec;
			for(unsigned i = 0; i < count; i++)
			{
				ss >> follower;
				tempVec.push_back(follower);
			}
			sort(tempVec.begin(), tempVec.end());

			intVectorMap[nodeId] = tempVec;
			// intVectorMap.push_back(tempVec);

			if(mapsize % 500000 == 0)
			{
				cout << "Reading File... Done with -- " << fileName << " -- " << mapsize << endl;
			}

			mapsize++;
		}

		intVectorFile.close();
		cout << "read -- " << fileName << endl;
	}
	else
	{
		cout << "Error opening file -- " << fileName << "\n";
	}
	return intVectorMap;
}

vector< vector <double> > readMultipleDoubleVectorsFromFile(string fileName)
{
	vector< vector <double> > vecVecdouble;
	ifstream vecVecdoubleFile;
	vecVecdoubleFile.open(fileName);

	string line;
	stringstream ss; 

	if(vecVecdoubleFile.is_open())
	{
		while(getline(vecVecdoubleFile, line))
		{
			ss.clear();
			ss.str("");

			ss << line;
			vector<double> tempVec;
			double tempVal;

			while(ss >> tempVal)
			{
				tempVec.push_back(tempVal);
			}

			vecVecdouble.push_back(tempVec);
		}

		vecVecdoubleFile.close();
	}
	else
	{
		cout << "Error opening file -- " << fileName << "\n";
	}
	return vecVecdouble;
}

vector <double>  readDoubleVectorFromFile(string fileName)
{
	vector <double> doubleVec;
	ifstream doubleVecFile;
	doubleVecFile.open(fileName);

	string line; 
   	stringstream ss; 

	ss.clear();
	ss.str("");

	if(doubleVecFile.is_open())
	{
		getline(doubleVecFile, line);

		ss << line;
		double tempVal;
		while(ss >> tempVal)
		{
			doubleVec.push_back(tempVal);
		}
	}
	else
	{
		cout << "Error opening file -- " << fileName << "\n";
	}

	return doubleVec;
}


map<ui, double> getUserBaseRates(string fileName)
{
	map<ui, double> ubRateMap;

	ifstream ubrFile;
	ubrFile.open(fileName);

	string line;
	stringstream ss;

	ui uid, twCount;			// time1, time2;
	double ubr;

	while(getline(ubrFile, line))
	{
		ss.clear();
		ss.str("");

		ss << line;

		ss >> uid >> twCount >> ubr;

		ubRateMap[uid] = ubr;
	}

	ubrFile.close();

	return ubRateMap;
}


int readHyperParametersVectorsFromFile()
{
	// hyperBeta = readDoubleVectorFromFile("hyperBetaValue.txt");
	// hyperGamma = readDoubleVectorFromFile("hyperGammaValue.txt");
	// hyperAlpha = readDoubleVectorFromFile("hyperAlphaValue.txt");


	// topicTopicProbVector = readMultipleDoubleVectorsFromFile("topicTopicProbVectors.txt");
	// userTopicPrefVector = readMultipleDoubleVectorsFromFile("nodeTopicProbVectors.txt");
	// topicWordProbVector = readMultipleDoubleVectorsFromFile("topicWordProbVectors.txt");

	vector<double> tempBeta(numTopics, 0.01);
	hyperBeta = tempBeta;

	vector<double> tempGamma(numTopics, 0.01);
	hyperGamma = tempGamma;

	vector<double> tempAlpha(vocabSize, alphaWord);
	hyperAlpha = tempAlpha;

	sumAlpha = accumulate(hyperAlpha.begin(), hyperAlpha.end(), 0.0);
	sumBeta = accumulate(hyperBeta.begin(), hyperBeta.end(), 0.0);
	sumGamma = accumulate(hyperGamma.begin(), hyperGamma.end(), 0.0);
	return 0;
}


int populateParentEventsForAllFromFile(string parentsFile, string parentExpFile)
{
	ifstream possParFile;
	possParFile.open(parentsFile);

	string line;
	stringstream ss;

	ui count, eid;
	if(possParFile.is_open())
	{

		while(getline(possParFile, line))
		{
			ss.clear();
			ss.str("");

			ss << line;

			ss >> count >> eid;

			vector<ui> parEvents;
			ui tempPar;
			for(ui i = 0; i < count; i++)
			{
				ss >> tempPar;
				parEvents.push_back(tempPar);
			}

			allPossibleParentEvents.push_back(parEvents);
			parEvents.clear();
		}

		possParFile.close();
	}
	else
	{
		cout << "Cannot open the possible parent events file\n";
		exit(0);
	}

	

	ifstream possExpFile;
	possExpFile.open(parentExpFile);

	if(possExpFile.is_open())
	{
		while(getline(possExpFile, line))
		{
			ss.clear();
			ss.str("");

			ss << line;

			ss >> count >> eid;
			vector <double> expEvents;
			double tempExp;
			for(ui i = 0; i < count; i++)
			{
				ss >> tempExp;
				expEvents.push_back(tempExp);
			}

			allPossibleParentEventsExponentials.push_back(expEvents);
			expEvents.clear();
		}

		possExpFile.close();
	}
	else
	{
		cout << "Cannot open the possible parent events Exponentials file\n";
		exit(0);
	}
	
	return 0;
}
