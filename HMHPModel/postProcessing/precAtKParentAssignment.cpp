/*
	Precision at k... for parent assignment...

	input file format:
	numOfProbValues eid  pid1 probval1 pid2 probval2 ...
	...
	...
	...
	...

*/

#include <iostream>
#include <cstdio>
// #include <boost/algorithm/string.hpp> 
// #include <boost/lexical_cast.hpp>
#include <fstream>
#include <map>
#include <unordered_map>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

using ui = unsigned int;


bool sortcol(const vector<double>& v1, const vector<double>& v2)
{
	// return v1[1] < v2[1];
	return v1[1] > v2[1]; 				// sort descending...
}


int main(int argc, char *argv[])
{	

	if(argc < 3)
	{
		cout << "Missing command line arguments...\n argv[1] = prec@kInputFile (Avg Probability values file...), argv[2] = prec@kOutputFile\n";
		exit(0);
	}

	// values of k
	// vector<int> kvalArr = {1, 3, 5, 10};
	vector<int> kvalArr = {0, 2, 4, 6, 9, 49, 99};

	ifstream origEventsFile;
	// origEventsFile.open("../centralFiles/allSyntheticEvents_Orig.txt");
	// origEventsFile.open("../events_semisyn_sample_gamma_1M.txt");
	origEventsFile.open("../inputFiles/eventsFile.txt");
	// allevents file...

	string line;
	stringstream ss;

	double evTime; 
	int evUid, evParId, evTopId, eventIndex, numOfProbValues;

	// map<int, int> evIdParIdMap;
	// unordered_map<int, int> evIdParIdMap;

	vector <int> evIdParIdMap;

	int linenum = 0;
	while(getline(origEventsFile, line))
	{
		ss.clear();
		ss.str("");

		ss << line;

		// should be changed...
		ss >> evTime >> evUid >> evParId >> evTopId;

		// evIdParIdMap[linenum] = evParId;
		evIdParIdMap.push_back(evParId);

		// cout << linenum << endl;
		if(linenum % 100000 == 0)
			cout << linenum << "\n";
		linenum++;
	}

	origEventsFile.close();

	cout << "Got original parent events\n";

	unordered_map <int, map<int, vector<bool> > >iterEidBoolVec;

	ifstream precAtKInputFile;
	// precAtKInputFile.open("precAtKInput.txt");
	precAtKInputFile.open(argv[1]);

	ofstream precAtKOutputFile;
	precAtKOutputFile.open(argv[2]);

	linenum = 0;

	// int numberOfIterVals;

	// int iterVal, numOfSamples;

	precAtKOutputFile << "lineNo. ";

	for(ui i = 0; i < kvalArr.size(); i++)
	{
		precAtKOutputFile << "#" << kvalArr[i] + 1 << " ";
	}

	precAtKOutputFile << "\n";

	string probLine;
	stringstream probStream;

	ofstream sortedValFile;
	sortedValFile.open("SortedParFile.txt");

	ofstream sortedValProbFile;
	sortedValProbFile.open("SortedProbFile.txt");

	while(getline(precAtKInputFile, probLine))
	{
		probStream.clear();
		probStream.str("");

		probStream << probLine;

		probStream >> numOfProbValues >> eventIndex; 

		int possParEid;
		double possParProbVal;

		vector< vector<double> > parIdProbVal;

		// while(numOfProbValues--)
		for(int j = 0; j < numOfProbValues; j++)
		{
			probStream >> possParEid >> possParProbVal;

			vector<double> tempVec;
			
			if(possParEid == linenum)
			{
				possParEid = -1;
			}

			tempVec.push_back(possParEid);
			tempVec.push_back(possParProbVal);

			parIdProbVal.push_back(tempVec);
			tempVec.clear();
		}

		sort(parIdProbVal.begin(), parIdProbVal.end(), sortcol);


		sortedValFile << parIdProbVal.size() << " " << linenum << " ";
		sortedValProbFile << parIdProbVal.size() << " " << linenum << " ";

		for(ui i = 0; i < parIdProbVal.size(); i++)
		{
			sortedValFile << parIdProbVal[i][0] << " ";
			sortedValProbFile << parIdProbVal[i][1] << " ";

		}

		sortedValFile << "\n";
		sortedValProbFile << "\n";

		// get the index of original parent... in predicted prob values...

		int origParent = evIdParIdMap[eventIndex];
		int origParIndex = 500010;

		for(ui l = 0; l < parIdProbVal.size(); l++)
		{
			if(parIdProbVal[l][0] == origParent)
			{
				origParIndex = l;
				break;
			}
		}

		precAtKOutputFile << linenum << " ";
		// various values of k for prec@k
		// vector<bool> presAbsOrigParentBoolVec;
		bool val;
		for(ui m = 0; m < kvalArr.size(); m++)
		{
			if(kvalArr[m] >= origParIndex)
			{
				val = 1;
			}
			else
			{
				val = 0;
			}
			// presAbsOrigParentBoolVec.push_back(val);
			precAtKOutputFile << val << " ";
		}
		precAtKOutputFile << "\n";

		// map<int, vector<bool> > eidBoolVec[linenum] = presAbsOrigParentBoolVec;

		// iterEidBoolVec =  eidBoolVec;

		linenum++;	
		if(linenum % 100000 == 0)
			cout << linenum << "\n";
	}

	sortedValFile.close();
	sortedValProbFile.close();

	precAtKOutputFile.close();
	precAtKInputFile.close();
	return 0;
}
