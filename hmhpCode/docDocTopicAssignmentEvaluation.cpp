/*
	Sampled document-document label switching...
*/

#include <iostream>
#include <cstdio>
#include <boost/algorithm/string.hpp> 
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <map>
#include <unordered_map>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

using ui = unsigned int;

vector<bool> calculateDocDocVector(vector<int> topicIdVec, string calling);

vector<bool> origDocDocVec;
vector<bool> predDocDocVec;

int main(int argc, char *argv[])
{
	if(argc < 2)
	{
		cout << "Missing Command line arguments...\nargv[1] = ModeTopicFile...\n";
		exit(0);
	}

	string line;
	stringstream ss;

	ifstream sampleSetFile;
	sampleSetFile.open("sampleSetIds.txt");

	map<int, bool> sampleSetMap;

	cout << "Get the sample of events ids..\n";
	int sampleId;

	while(getline(sampleSetFile, line))
	{
		ss.clear();
		ss.str("");

		ss << line;

		ss >> sampleId;

		sampleSetMap[sampleId] = 1;
	}

	cout << "Sample set size = " << sampleSetMap.size() << "\n";
	sampleSetFile.close();

	ifstream origEventsFile;
	// origEventsFile.open("../centralFiles/allSyntheticEvents_Orig.txt");
	// origEventsFile.open("../centralFiles/allSemiSyntheticEvents_global_C5_AlphaExp_1.0.txt");
	origEventsFile.open("../inputFiles/eventsFile.txt");
	// allevents file...

	double evTime; 
	int evUid, evParId, evTopId, eventIndex;

	vector<int> origTopicIdVec;

	int linenum = 0;

	map<int, int> origDistinctTopicMap;


	if(origEventsFile.is_open())
	{
		while(getline(origEventsFile, line))
		{	
			if(sampleSetMap.find(linenum) != sampleSetMap.end())
			{
				// cout << line << endl;
				ss.clear();
				ss.str("");

				ss << line;

				ss >> evTime >> evUid >> evParId >> evTopId;

				origTopicIdVec.push_back(evTopId);

				origDistinctTopicMap[evTopId] = 1;
			}
			linenum++;

		}
		origEventsFile.close();

		cout << "Got the topic Ids for the original set of events...\n";
		
	}
	else
	{
		cout << "cannot open original events file...\n";
		exit(0);
	}

	cout << "Number of distinct events = " << origTopicIdVec.size() << "\n";
	cout << "original distinct topics = " << origDistinctTopicMap.size() << "\n";

	ifstream topicAssignFile;
	topicAssignFile.open(argv[1]);
	// allevents file...

	vector<int> predTopicIdVec;

	linenum = 0;
	int predTopic;

	map<int, int> predDistinctTopicMap;

	while(getline(topicAssignFile, line))
	{
		if(sampleSetMap.find(linenum) != sampleSetMap.end())
		{
			ss.clear();
			ss.str("");

			ss << line;

			ss >> eventIndex >> predTopic;

			predTopicIdVec.push_back(predTopic);

			predDistinctTopicMap[predTopic] = 1;
		}

		linenum++;
	}
	topicAssignFile.close();

	cout << "Number of distinct events = " << predTopicIdVec.size() << "\n";
	cout << "Predicted distinct topics = " << predDistinctTopicMap.size() << "\n";

	cout << "Got the topic Ids for the predicted set of events...\n";

	origDocDocVec = calculateDocDocVector(origTopicIdVec, "Original Set");
	predDocDocVec = calculateDocDocVector(predTopicIdVec, "Predicted Set");

	double count = 0;

	int TP = 0, TN = 0, FP = 0, FN = 0;

	for(unsigned int i = 0; i < origDocDocVec.size(); i++)
	{
		if(origDocDocVec[i] == predDocDocVec[i])
		{
			count++;
		}

		if(origDocDocVec[i] == 1 && predDocDocVec[i] == 1)
		{
			TP += 1;
		}
		else if(origDocDocVec[i] == 1 && predDocDocVec[i] == 0)
		{
			FN += 1; 
		}
		else if(origDocDocVec[i] == 0 && predDocDocVec[i] == 1)
		{
			FP += 1;
		}
		else
		{
			TN += 1;
		}
	}


    double prec = (TP * 1.0 )/ (TP + FP);
    double recall = (TP * 1.0 )/ (TP + FN);
    double f1score = 2 * (prec * recall) / (prec + recall);

    cout << "Doc Doc Comparison match = " << (count * 1.0)/ origDocDocVec.size() << "\n";
    cout << "Precision = " << prec << endl;
    cout << "Recall = " << recall << endl;
    cout << "F1 score = " << f1score << endl;

    cout << "FP = " << FP << endl;
    cout << "FN = " << FN << endl;
    return 0;
}


vector<bool> calculateDocDocVector(vector<int> topicIdVec, string calling)
{ 
	vector<bool> docDocVec;

	int zeroCount = 0, oneCount = 0;

	for(unsigned int i = 0; i < topicIdVec.size(); i++)
	{
		int countj = 0;
		for(unsigned int j = i+1; j < topicIdVec.size(); j++)
		{
			if(topicIdVec[i] == topicIdVec[j])
			{
				docDocVec.push_back(1);
				oneCount++;
			}
			else
			{
				docDocVec.push_back(0);	
				zeroCount++;
			}
			countj++;
		}
		// cout << countj << " ";
	}
	// cout << "\n";
	cout << calling << " zeroCount = " << zeroCount << " oneCount = " << oneCount << endl;

	return docDocVec;
}
