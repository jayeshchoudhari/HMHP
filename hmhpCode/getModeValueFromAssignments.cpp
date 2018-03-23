/*
	Get predicted mode value (occurring highest number of times) for each event from the assignment file with iterations...
	Can be used for both parent assignment as well as topic assignment...
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
	if(argc > 2)
	{
		ifstream assignmentFile;
		assignmentFile.open(argv[1]);

		map<ui, int> eidAssignmentMap;

		ui eid;

		string line;
		stringstream ss;

		int linenum = 0;
		while(getline(assignmentFile, line))
		{
			ss.clear();
			ss.str("");

			ss << line;

			ss >> eid;

			int assignment;
			unordered_map<int, int> assignmentCountMap;
			int maxCount = 0;
			int maxCountAssign = -1;

			while(ss >> assignment)
			{
				assignmentCountMap[assignment]++;
				if(assignmentCountMap[assignment] > maxCount)
				{
					maxCount = assignmentCountMap[assignment];
					maxCountAssign = assignment;
				}
			}

			eidAssignmentMap[linenum] = maxCountAssign;

			linenum++;
		}

		assignmentFile.close();

		// write mode assignments to the file -- eid assignment

		ofstream outAssignFile;
		outAssignFile.open(argv[2]);

		map<ui, int>::iterator eidAssignmentMapIterator;

		for(eidAssignmentMapIterator = eidAssignmentMap.begin(); eidAssignmentMapIterator != eidAssignmentMap.end(); eidAssignmentMapIterator++)
		{
			outAssignFile << eidAssignmentMapIterator->first << " " << eidAssignmentMapIterator->second << "\n"; 
		}

		outAssignFile.close();
	}
	else
	{
		cout << "Missin command line arguments -- argv[1] should be assignmentFile with 10 or more values and argv[2] must be the output file..\n";  
	}
	

	return 0;
}
