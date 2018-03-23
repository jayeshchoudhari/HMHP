/*
	Get top-50 possible/probable parent events for all the events...
	And store them in a file... so that we do not calculate them next time, but just populate from the file...
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
using li = long int;

vector < vector <ui> > reverseFollowersMap;
vector< vector <li> > newSyntheticEvents;
unordered_map<ui, double> eventIndexTimestamps;

vector <vector <ui> > allPossibleParentEvents;
vector <vector <double> > allPossibleParentEventsExponentials;
vector <vector <double> > allPossibleParentTimeDiff;
unordered_map<ui, vector<ui> > nodeEventsMap;

ui maxNumNodes = 1000002;					// 7697889	

// double timeKernelMultiplier = 1;
// double alphaMultiplier = 0.1
double timeKernelMultiplier = 1;

vector < vector <ui> > readIntVectorMapFromFile(string fileName);
vector< vector <li> > getSyntheticEventsFromFile(string fileName);

ui getFirstIndexGreaterThanEventIndex(vector <ui> neighborNodeEvents, ui eventIndex);
int populateParentEventsForAll();

bool sortcol(const vector<double>& v1, const vector<double>& v2)
{
	// return v1[1] < v2[1];
	return v1[1] > v2[1]; 				// sort descending...
}

ofstream allPossParentEventsFile, allPossExpFile, allPossTimeDiffFile;

int maxAllowedTime = 13;
ui maxCandParents = 100;

int main(int argc, char *argv[])
{
	reverseFollowersMap = readIntVectorMapFromFile("../mapped_users_friends_restricted_tweets_top_5k_hashes.txt");
	// reverseFollowersMap = readIntVectorMapFromFile("../friends_BA_100_40.0.txt");
	cout << "Got the reverse followers vector... " << reverseFollowersMap.size() << "\n";

    //newSyntheticEvents =  getSyntheticEventsFromFile("events_100exp_0.1.txt");
    newSyntheticEvents =  getSyntheticEventsFromFile(argv[1]);
	cout << "Got all the events -- " << newSyntheticEvents.size() << "\n";
	// newSyntheticEvents =  getSyntheticEventsFromFile("oscars_tweets_events.txt");

	//allPossParentEventsFile.open("allParentsForHTEvents_0.10_timediff_10_exact_cand_parents.txt");
	allPossParentEventsFile.open(argv[2]);
	allPossExpFile.open(argv[3]);

	cout << "Getting top 100 possible parent events for all the events...with time diff upto 200\n";
	populateParentEventsForAll();

	cout << "Writing parent events and exponentials to the file...\n";


	allPossParentEventsFile.close();
	allPossExpFile.close();
	// allPossTimeDiffFile.close();
	
	return 0;
}

int populateParentEventsForAll()
{
	// These neighbors are to be found on the reverseMap...

	// get indices of all the possible parent events...
	// parent event can come from neighbors only...
	// get the list of the neighbors

	double avgPossParents = 0;

	vector < vector <ui> > tempParentEventForAll;

	double countAvgParents = 0;

	// #pragma omp parallel for 
	// #pragma omp parallel num_threads(8)
	for(ui i = 0; i < newSyntheticEvents.size(); i++)
	// for(ui i = 0; i < 20; i++)
	{
		ui eventIndex = i;
		ui eventNode = newSyntheticEvents[i][1];
        // int eventParent = newSyntheticEvents[i][2];
		double eventTime = eventIndexTimestamps[i];

		vector <ui> neighbors = reverseFollowersMap[eventNode];
		
		vector < vector <double> > parentExpValues;

		vector <vector <double> > parentDiffValues;
        
        // if(eventParent == eventNode || eventParent < 0)
                // neighbors.clear();

		// vector <ui> tempVec;

		if (neighbors.size() > 0)
		{
			for(unsigned int k = 0; k < neighbors.size(); k++)
			{
				ui neighborNode = neighbors[k];
				// get its events of lesser timestamp
				vector<ui> neighborNodeEvents = nodeEventsMap[neighborNode];
				if(neighborNodeEvents.size() > 0)
				{
					// cout << neighborNode << endl;
					int startPoint;

					if(neighborNodeEvents[neighborNodeEvents.size() - 1] < eventIndex)
					{
						startPoint = neighborNodeEvents.size() - 1;
					}
					else if(neighborNodeEvents[0] > eventIndex)
					{
						startPoint = -1;
					}
					else
					{
						startPoint = getFirstIndexGreaterThanEventIndex(neighborNodeEvents, eventIndex);
					}

					// short int count_2 = 0;
					// for(ui j = 0; j < neighborNodeEvents.size(); j++)
					for(int j = startPoint; j >= 0; j--)
					{	
						vector<double> tempVec;
						//if(neighborNodeEvents[j] < eventIndex && neighborNodeEvents[j] >= eventParent)
                        if(neighborNodeEvents[j] < eventIndex)
						{
							double thisEventTime = eventIndexTimestamps[neighborNodeEvents[j]];
							// possibleParentEvents.push_back(neighborNodeEvents[j]);
							double timeDiff = (eventTime - thisEventTime);

							if(timeDiff < maxAllowedTime)
							{	
								double expVal = exp(-timeKernelMultiplier * timeDiff);
								// possibleParentExp.push_back(expVal);

								tempVec.push_back(neighborNodeEvents[j]);
								tempVec.push_back(expVal);

								// tempVec.push_back(eventTime - thisEventTime);

								parentExpValues.push_back(tempVec);
								tempVec.clear();
							}
							else
							{
								break;
							}
						}
					}
				}
			}
		}

		// avgPossParents += possibleParentEvents.size();
		
		vector <ui> possibleParentEvents; 
		vector <double> possibleParentExp, possParentTimeDiff;
		
		if(parentExpValues.size() > 0)
		{
			sort(parentExpValues.begin(), parentExpValues.end(), sortcol);

			ui getTopEvents = parentExpValues.size() > maxCandParents ? maxCandParents : parentExpValues.size();
			//ui getTopEvents = parentExpValues.size();

			countAvgParents += getTopEvents;

			allPossParentEventsFile << getTopEvents << " " << i << " ";
			allPossExpFile << getTopEvents << " " << i << " ";

			for(ui j = 0; j < getTopEvents; j++)
			{
				// possibleParentEvents.push_back((ui)parentExpValues[j][0]);
				// possibleParentExp.push_back(parentExpValues[j][1]);
				// possParentTimeDiff.push_back(parentExpValues[j][2]);
				
				allPossParentEventsFile << parentExpValues[j][0] << " ";
				allPossExpFile << parentExpValues[j][1] << " ";
			}
		}
		else
		{
			allPossParentEventsFile << parentExpValues.size() << " " << i;
			allPossExpFile << parentExpValues.size() << " " << i;			
		}

		allPossParentEventsFile << "\n";
		allPossExpFile << "\n";

		// tempParentEventForAll.push_back(possibleParentEvents);
		// allPossibleParentEvents.push_back(possibleParentEvents);
		// allPossibleParentEventsExponentials.push_back(possibleParentExp);
		// allPossibleParentTimeDiff.push_back(possParentTimeDiff);

		// avgPossParents += possibleParentEvents.size();

		// possibleParentExp.clear();
		// possibleParentEvents.clear();

		if(i%100000 == 0)
			cout << "Got parents for " << i << "  events \n";
	}

	// cout << "Average Possible Parents -- " << avgPossParents / newSyntheticEvents.size() << endl;
	cout << "Average Possible Parents -- " << countAvgParents / newSyntheticEvents.size() << endl;

	// return tempParentEventForAll;
	return 0;
}

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

			// if(mapsize == 500000)
			// {
			// 	intVectorFile.close();			
			// 	break;
			// }
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

ui getFirstIndexGreaterThanEventIndex(vector <ui> neighborNodeEvents, ui eventIndex)
{

	ui low = 0, high = neighborNodeEvents.size();
	ui mid;

	while(low != high)
	{
		mid = (low + high) / 2;

		if(neighborNodeEvents[mid] < eventIndex)
		{
			low = mid + 1;
		}
		else
		{
			high = mid;
		}
	}

	return low;
}


vector< vector <li> > getSyntheticEventsFromFile(string fileName)
{
	// vector < vector <double> > allEvents;
	vector < vector <li> > allEvents;

	ifstream allEventsFile;
    allEventsFile.open(fileName);

    string line;
    stringstream ss;

    unsigned i = 0;

    if(allEventsFile.is_open())
    {
    	while(getline(allEventsFile, line))
		{
			double eventTime;
			li eventNode, eventParent, eventTopic;

			ss.clear();
			ss.str("");

			ss << line;
			ss >> eventTime >> eventNode >> eventParent >> eventTopic;
			// ss >> eventTime >> eventNode >> eventParent;

			eventIndexTimestamps[i] = eventTime;

			vector<li> tempEvent;
			tempEvent.push_back(i);						// eventIndex...
			tempEvent.push_back(eventNode);
            tempEvent.push_back(eventParent);

			allEvents.push_back(tempEvent);

			nodeEventsMap[eventNode].push_back(i);

			line.clear();

			i++;

			tempEvent.clear();
		}
    }
	else
	{
		cout << "Error opening file -- " << fileName << "\n";
	}
	cout << "read " << allEvents.size() << " events\n";

	cout << "First event time -- " << eventIndexTimestamps[0] << "\n";

	allEventsFile.close();
	return allEvents;
}
