/*
	Calculate Predicted Accuracy...

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


int main(int argc, char *argv[])
{
	if(argc < 2)
	{
		cout << "Missing command line argument.. argv[1] = parentModeAssignmentFile\n";
		exit(0);
	}

	string line;
	stringstream ss;

	/*
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

	sampleSetFile.close();
	
	cout << "Got the set of sample events...\n";
	*/
	ifstream origEventsFile;
	// origEventsFile.open("../centralFiles/allSyntheticEvents_Orig.txt");
	// origEventsFile.open("../events_semisyn_sample_gamma_1M.txt");
	origEventsFile.open("../inputFiles/eventsFile.txt");
	// allevents file...

	double evTime; 
	int evUid, evParId, evTopId, eventIndex;

	map<int, int> evIdParIdMap;

	int linenum = 0;

	while(getline(origEventsFile, line))
	{
		// if(sampleSetMap.find(linenum) != sampleSetMap.end())
		// {
			ss.clear();
			ss.str("");

			ss << line;

			ss >> evTime >> evUid >> evParId >> evTopId;

			evIdParIdMap[linenum] = evParId;
		// }

		linenum++;
	}
	origEventsFile.close();

	cout << "Got the parent Ids for the original set of events...\n";

	ifstream parAssignFile;
	// parAssignFile.open("parentAssignmentHighestPercent.txt");
	parAssignFile.open(argv[1]);

	linenum = 0;

	int assignedParId, origParId;
	double count = 0;

	ofstream unmatchedModeFile;
	unmatchedModeFile.open("unmatchedMode.txt");


	while(getline(parAssignFile, line))
	{
		// if(sampleSetMap.find(linenum) != sampleSetMap.end())
		// {
			ss.clear();
			ss.str("");

			ss << line;

			ss >> eventIndex >> assignedParId;

			origParId = evIdParIdMap[linenum];
					
			// cout << origParId << " " << assignedParId << "\n";

			if(origParId == assignedParId)
			{
				count += 1;
			}
			else
			{
				unmatchedModeFile << linenum << endl;
			}
		// }

		linenum++;
	}

	unmatchedModeFile.close();
	// double predictedAcc = (count*1.0) / sampleSetMap.size();
	double predictedAcc = (count*1.0) / linenum;

	cout << "Predicted Accuracy = " << predictedAcc << endl;

	return 0;
}
