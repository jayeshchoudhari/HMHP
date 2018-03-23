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
#include <algorithm>
#include <numeric>	
#include <random>
#include <chrono>
#include <cfenv>


// #pragma STDC FENV_ACCESS ON

using namespace std;
using namespace std::chrono;

// using ll = unsigned long long int;
using ll = long long int;
using li = long int;
using ui = unsigned int;

#define INF std::numeric_limits<int>::max()					//use limits for infinity
#define PI 3.14159265

// function declaration

// Read data from files....
vector < vector <li> > getSyntheticDocsFromFile(string fileName);
int readHyperParametersVectorsFromFile();

// initialization
vector< vector <ll> > initializeForSampler(vector< vector <double> > allEvents, int flag);
int initializeAvgTopicProbabilityVectors();

// sample topic
int sampleTopicAssignment(int ITE);
// ll getSampledTopicAssignment(vector<li> doc, ll eid);
ll getSampledTopicAssignment(vector<li> doc, ll eid, int ITE);
double getThirdTermOfTopicAssignmentCondProb(ll topic, vector<li> doc, vector<unsigned long> wordHist);


int createWordHistForAllDocs();
unordered_map <li, li> getHistOfWordsOverWordsFromDoc(vector<li> doc);

// decreament
int decreamentCountFromMatrices(ll eventTopic, vector <li> doc);


// increament
int increamentCountToMatrices(ll eventTopic, vector<li> doc);


// validations
// int countCorrectFractionTopicAssignments();
int writeEvery10ItersTopicAssignmentToFile();
int writeAvgTopicProbVectorsToFile();
int writeTopicWordDistributionToFile();


// util functions
int getSampleFromDiscreteDist(vector<double> normalizedProbVector);

// global members declaration

vector< vector <li> > allDocs;
vector <ll> assignedSampleTopics;

double sumAlpha;

// Count Matrices... 
vector < vector <unsigned long> > wordHistAllDocsVector;

// validation Maps
map <int, vector<int> > topicDistEvery10thIter;

int BURN_IN = 200;
ui ITERATIONS = 301;

ui numTopics = 100;

int vocabSize = 33900; 									// 114462
// int vocabSize = 28500; 								// 28101
unsigned long totalWords = 0;

double alphaWord = 0.1;
double betaTopic = 0.1;

//std::vector<std::vector<int>> vec_2d(rows, std::vector<int>(cols, 0));
vector < vector < unsigned long > > NTWCountVec(numTopics, vector < unsigned long >(vocabSize, 0));
vector < unsigned long > NTWSumWordsVec(numTopics);

// vector < vector <double> > avgTopicProbVector;

vector<int> topicPopularityCount (numTopics, 0);

int main(int argc, char *argv[])
{
    if(argc > 1)
	{
		BURN_IN = atoi(argv[1]);
		ITERATIONS = atoi(argv[2]);
	}

	cout << "Will be running sampler for BURN_IN = " << BURN_IN << " and ITERATIONS = " << ITERATIONS << "\n";

	// newSyntheticEvents = getSyntheticEventsFromFile("semiSyntheticTrainingEvents.txt");
	allDocs = getSyntheticDocsFromFile("../docs_semisyn_sample_gamma_1M.txt");
	cout << "Got all the docs..." << allDocs.size() << "\n";
	// allDocs = getSyntheticDocsFromFile("test-semi-docs.txt");
	
	// initializeAvgTopicProbabilityVectors();

	cout << "Creating hist of words for each doc...\n";
	createWordHistForAllDocs();
	cout << "Created hist of words for each doc...\n";

	for(unsigned int i = 0; i < allDocs.size(); i++)
	{
		// assignedSampleTopics.push_back(-1);
		int initTopic = i % numTopics;
		assignedSampleTopics.push_back(initTopic);
		topicPopularityCount[initTopic]++;

		for(ui j = 0; j < allDocs[i].size(); j++)
		{
			int word = allDocs[i][j];
			NTWCountVec[initTopic][word]++;
		}

		NTWSumWordsVec[initTopic] += allDocs[i].size();

	}

	sumAlpha = alphaWord * vocabSize;
	if(sumAlpha < 0)
	{
		cout << "wrong vocabSize \n";
		exit(0);
	}

	cout << "Sampling Topics for Documents\n";

	for(unsigned int ITE = 0; ITE < ITERATIONS; ITE++)
	{
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		sampleTopicAssignment(ITE);
    	high_resolution_clock::time_point t2 = high_resolution_clock::now();
		
		auto duration = duration_cast<microseconds>( t2 - t1 ).count();

		cout << "ITE -- " << ITE << " Duration -- " << duration << endl;

	}	// ITERATIONS...

	writeEvery10ItersTopicAssignmentToFile();
	// writeAvgTopicProbVectorsToFile();
	writeTopicWordDistributionToFile();

	return 0;
}

///////////////////////////////
///////// INITIALIZATION //////
///////////////////////////////

/*
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
*/

int createWordHistForAllDocs()
{
	for(unsigned int i = 0; i < allDocs.size(); i++)
	{
		vector<long int> doc = allDocs[i];

		unordered_map<li, li> wordHist;
		unordered_map<li, li>::iterator wordHistIt;
		wordHist = getHistOfWordsOverWordsFromDoc(doc);

		// wordHistAllDocs[i] = wordHist;
		vector <unsigned long> tempWordCounts;
		for(wordHistIt = wordHist.begin(); wordHistIt != wordHist.end(); wordHistIt++)
		{
			tempWordCounts.push_back(wordHistIt->first);
			tempWordCounts.push_back(wordHistIt->second);
		}
		wordHistAllDocsVector.push_back(tempWordCounts);
	}

	cout << "Got the Hist of words for all Docs\n";

	return 0;
}


///////////////////////////////
///////     OUTPUTS    	 //////
///////////////////////////////

int writeEvery10ItersTopicAssignmentToFile()
{
	ofstream topicAssignmentsFile;
	topicAssignmentsFile.open("B1_topicAssignments_LDAMM_100_Topics_300_ites.txt");

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


/*
int writeAvgTopicProbVectorsToFile()
{
	ofstream avgTopicProbVecFile;
	// avgProbVecFile.open("avgParProb_OurModel.txt");
	avgTopicProbVecFile.open("B1_avgTopicProb_LDAMM.txt");

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
*/

int writeTopicWordDistributionToFile()
{
	ofstream tWDistFile;
	// tWDistFile.open("topicWordDistFile_1M.txt");
	tWDistFile.open("topicWordCount_1M_50_Topics.txt");

	for(unsigned int i = 0; i < NTWCountVec.size(); i++)
	{
		int sumCount = NTWSumWordsVec[i];
		for(unsigned int j = 0; j < NTWCountVec[i].size(); j++)
		{
			int count = NTWCountVec[i][j];
			if(sumCount > 0)
			{	
				tWDistFile << count << " ";
			}
			else
			{
				tWDistFile << 0 << " ";	
			}
		}	

		tWDistFile << "\n";
	}

	tWDistFile.close();

	return 0;
}

///////////////////////////////
///////  SAMPLE TOPIC    //////
///////////////////////////////

int sampleTopicAssignment(int ITE)
{
	// cout << "Sample Topic Assignment\n";
	for(unsigned int i = 0; i < allDocs.size(); i++)
	// for(unsigned int i = 0; i < 275; i++)
	{
		decreamentCountFromMatrices(assignedSampleTopics[i], allDocs[i]);
		ll assignedTopic =  getSampledTopicAssignment(allDocs[i], i, ITE);
		increamentCountToMatrices(assignedTopic, allDocs[i]);
		assignedSampleTopics[i] = assignedTopic;

		if(ITE > BURN_IN && ITE % 10 == 0)
		{
			topicDistEvery10thIter[i].push_back(assignedTopic);
		}
	}

	return 0;
}


ll getSampledTopicAssignment(vector<li> doc, ll eid, int ITE)
{
	ll assignedTopic = -1;

	vector<unsigned long> wordHistVec = wordHistAllDocsVector[eid];
	vector <double> calculatedProbVec (numTopics, 0.0);

	if(wordHistVec.size() > 0)
	{
		double normalization = 0.0;

		for(ll topic = 0; topic < numTopics; topic++)
		{
			double thirdTerm = 0, firstTerm = 0;
			firstTerm = betaTopic + topicPopularityCount[topic];
			thirdTerm = getThirdTermOfTopicAssignmentCondProb(topic, doc, wordHistVec);
			calculatedProbVec[topic] = firstTerm * thirdTerm;
			normalization += (firstTerm * thirdTerm);
		}

		double normalizationFactor = 1.0/normalization;
		for( unsigned i = 0; i < calculatedProbVec.size(); i++)
		{
			calculatedProbVec[i] *= normalizationFactor; 
		}
		assignedTopic = getSampleFromDiscreteDist(calculatedProbVec);

		/*
		// updating for avgProbability of the topic for each event... at the end lets write it to a file...  
		if(ITE > BURN_IN && ITE%10 == 0 && assignedTopic > -1)
		{
			for(ui i = 0; i < avgTopicProbVector[eid].size(); i++)
			{
				avgTopicProbVector[eid][i] += calculatedProbVec[i];
			}
		}
		*/
	}
	else
	{
		assignedTopic = eid % numTopics;
	}
	calculatedProbVec.clear();

	return assignedTopic;
}


double getThirdTermOfTopicAssignmentCondProb(ll topic, vector<li> doc, vector<unsigned long> wordHistVec)
{

	long double finalThirdValue = 0;
	long double finalThirdValueNume = 1;
	long double finalThirdValueDenom = 1;

	unsigned i = 0;

	while(i < wordHistVec.size())
	{
		unsigned long word = wordHistVec[i];
		i++;
		unsigned long wordCount = wordHistVec[i];

		// as for all the words this is same... we do not really need hyperAlpha right now...
		double baseterm = alphaWord + NTWCountVec[topic][word];

		long double valueForWord = 1.0;
		// double valueForWord = 1.0;

		for(unsigned int i = 0; i < wordCount; i++)
		{
			valueForWord *= (baseterm + i);
		}

		finalThirdValueNume *= valueForWord;
		i++;
	}

	double baseTermDenom = sumAlpha + NTWSumWordsVec[topic];

	for(unsigned int i = 0; i < doc.size(); i++)
	{
		finalThirdValueDenom *= (baseTermDenom + i);
	}

	finalThirdValue = finalThirdValueNume / finalThirdValueDenom;

	return finalThirdValue;
}


// this is required for calculating the lastterm of the probability for topic assignment...
unordered_map <li, li> getHistOfWordsOverWordsFromDoc(vector<li> doc)
{
	unordered_map <li, li> wordCounts;

	for(unsigned int i = 0; i < doc.size(); i++)
	{
		// ll word = doc[i];
		wordCounts[doc[i]] += 1;
	}

	return wordCounts;
}

// Decreament Counters
int decreamentCountFromMatrices(ll eventTopic, vector <li> doc)
{
	// cout << "decreament count from topic-topic matrix.. from the cell parentTopic -> eventTopic";
	// decreament only if the eventTopic is valid...
	if(eventTopic >= 0)
	{
		if(topicPopularityCount[eventTopic] > 0)
		{
			topicPopularityCount[eventTopic] -= 1;
		}
		else
		{
			cout << "Some Issue in counting popularity of topics...\n";
		}

		if(NTWSumWordsVec[eventTopic] >= doc.size())
		{
			NTWSumWordsVec[eventTopic] -= doc.size();
		}
		else
		{	
			cout << "Somethings wrong with the Topic Word Counts... \n";
			NTWSumWordsVec[eventTopic] = 0;	
		}

		for(unsigned int i = 0; i < doc.size(); i++)
		{
			if(NTWCountVec[eventTopic][doc[i]] > 0)
			{
				NTWCountVec[eventTopic][doc[i]] -= 1;
			}
		}
	}
	return 0;
}

// Increament Counters
int increamentCountToMatrices(ll eventTopic, vector<li> doc)
{
	// increament counts in topic-word...
	// increamentCountInTopicWord(eventTopic, doc);

	NTWSumWordsVec[eventTopic] += doc.size();

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

	topicPopularityCount[eventTopic] += 1;

	return 0;
}

///////////////////////////////
///////  UTIL FUNCTIONS  //////
///////////////////////////////

int getSampleFromDiscreteDist(vector<double> normalizedProbVector)
{
	random_device rd;
	mt19937 gen(rd());

	discrete_distribution<int> distribution(normalizedProbVector.begin(), normalizedProbVector.end());
	int ind = distribution(gen);

	return ind;
}

///////////////////////////////
///////// READING DATA ////////
//////// AND PARAMETERS ///////
//////// FROM FILES ///////////
///////////////////////////////

vector < vector <li> > getSyntheticDocsFromFile(string fileName)
{
	vector < vector <li> > allDocs;

	double avgDocLength = 0;

	ifstream allDocsFile;
    allDocsFile.open(fileName);

    string line;
    stringstream ss;
    ll i = 0;
    while(getline(allDocsFile, line))
	{
		ss.clear();
		ss.str("");

		ss << line;

		int word;
		vector <li> currDoc;
		while(ss >> word)
		{
			currDoc.push_back(word);
		}
		
		// eventIndexDocSizeMap[i] = currDoc.size();
		i++;

		sort(currDoc.begin(), currDoc.end());

		avgDocLength += currDoc.size();

		allDocs.push_back(currDoc);
		line.clear();
	}

	totalWords = avgDocLength;
	avgDocLength = avgDocLength/allDocs.size();	

	cout << "Total Number of words = " << totalWords << "\n";
	cout << "Read " << allDocs.size() << " documents\n";
	cout << "Avg Doc Length = " << avgDocLength << endl;

	allDocsFile.close();
	return allDocs;
}
