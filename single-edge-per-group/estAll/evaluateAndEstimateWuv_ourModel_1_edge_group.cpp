// estimate Wuv per edge in the predicted set of parents... 

#include <iostream>				//for basic C++ functions 
#include <cstdio>				//for std C functions
#include <fstream>				//for i/o stream
#include <sstream>				//for parse data using stringstreams...
#include <string>				//for C++ string functions
#include <cstring>				//we need this for memset... and string functions from C
#include <unordered_map>		//to maintain the map of user-tweet_count, and user_user_tweet_count
#include <map>					//to maintain the map of user-tweet_count, and user_user_tweet_count
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

using ui = unsigned int;

unordered_map <ui, unordered_map <ui, double> > readUserUserInfluenceFromFile(string fileName);
double calcVar(vector<double> errorVec, double meanVal);

unordered_map <ui, unordered_map<ui, double> > userUserInfluence;
unordered_map <ui, unordered_map<ui, double> > predUserUserInfluence;
unordered_map <ui, unordered_map<ui, int> > edgeMarkWuv;
unordered_map <ui, unordered_map<ui, int> > predNuv;

double alphaDash = 0.01, betaDash = 1;

int main()
{
	ifstream eventsFile;
	eventsFile.open("../events_semisyn_sample_gamma_1M.txt");

	string line;
	stringstream ss;

	double evTime;
	int evUid, evParId, evTopic;

	vector <int> parentNodeVec;

	unordered_map <ui, ui> NuValuesVec;

	int linenum = 0;
	int countEdges = 0;
	while(getline(eventsFile, line))
	{
		ss.clear();
		ss.str("");

		ss << line;

		ss >> evTime >> evUid >> evParId >> evTopic;
		parentNodeVec.push_back(evUid);

		// if(linenum < 150000)
		// {
		NuValuesVec[evUid] += 1;
		// }

		if(evParId > -1)
		{
			int parNode = parentNodeVec[evParId];
			edgeMarkWuv[parNode][evUid] += 1;
			countEdges++;
		}

		linenum++;
	}

	eventsFile.close();

	cout << "Marked required edges..." << edgeMarkWuv.size() << " Total Edges = " << countEdges << "\n";

	cout << "Getting user user influence...\n";
	userUserInfluence = readUserUserInfluenceFromFile("../../mapped_users_influence_restricted_tweets_top_5k_hashes_Gamma_alpha_beta_Wuv_0.01.txt");

	cout << "Got the user user influence map... " << userUserInfluence.size() << "\n";

	// cout << "Getting predicted user user influence...\n";
	// predUserUserInfluence = readUserUserInfluenceFromFile("./B5_userUserInf_OurModel_estAll.txt");
	// cout << "Got the predicted user user influence map... " << predUserUserInfluence.size() << "\n";

	// get the predicted wuv using (Nuv+alphaDash) / (Nu+betaDash)

	int pid, evId;

	// vector<int> predParentIdsVec;
	ui childNode;

	ifstream predParentsFile;
	predParentsFile.open("./lastIterationParents.txt");

	linenum = 0;
	if(predParentsFile.is_open())
	{
		while(getline(predParentsFile, line))
		{
			ss.clear();
			ss.str("");

			ss << line; 
			ss >> evId >> pid;

			// predParentIdsVec.push_back(pid);

			if(pid > -1)
			{
				ui parNode = parentNodeVec[pid];
				childNode = parentNodeVec[linenum];

				predNuv[parNode][childNode] += 1;
			}

			linenum += 1;
		}
	}
	else
	{
		cout << "cannot open pred parents file \n";
		exit(0);
	}

	cout << "Got the predicted Nuv values...\n";
	cout << "Calculating the Wuv values...\n";

	unordered_map <ui, unordered_map<ui, int> >::iterator edgeMarkWuvIt;

	ofstream predValuesFile;
	predValuesFile.open("origPredWuvValues.txt");

	// predValuesFile << "origUVInf" << " " << "predUVInf" << " " << "Transactions" << " " << "currSOSError" << " " << "currAbsError" << " " << "currMapeError" << "\n";
	predValuesFile << "origUVInf" << " " << "predUVInf" << " " << "Nuv" << " " << "Nu" << " " << "ape" << "\n";

	for(edgeMarkWuvIt = edgeMarkWuv.begin(); edgeMarkWuvIt != edgeMarkWuv.end(); edgeMarkWuvIt++)
	{
		int uNode = edgeMarkWuvIt->first;
		unordered_map<ui, int> interMap = edgeMarkWuvIt->second;

		unordered_map<ui, int>::iterator interMapIt;

		for(interMapIt = interMap.begin(); interMapIt != interMap.end(); interMapIt++)
		{
			int vNode = interMapIt->first;

			double origUVInf = userUserInfluence[uNode][vNode];

			double predUVInf = ((predNuv[uNode][vNode] + alphaDash) * 1.0 ) / (NuValuesVec[uNode] + betaDash);

			double ape = ((abs(origUVInf - predUVInf)) * 1.0) / origUVInf;

			predValuesFile << origUVInf << " " << predUVInf << " " << predNuv[uNode][vNode] << " " << NuValuesVec[uNode] << " " << ape << "\n";
			
		}
	}

	predValuesFile.close();

	
	return 0;

}


double calcVar(vector<double> errorVec, double meanVal)
{
	double varVal = 0;

	double sumVal = 0.0;
    double temp = 0.0;
   
    for (unsigned int j =0; j < errorVec.size(); j++)
    {
        temp = pow((errorVec[j] - meanVal),2);
        sumVal += temp;
    }
   
    varVal = sumVal/(errorVec.size() - 1);

	return varVal;
}


unordered_map <ui, unordered_map <ui, double> > readUserUserInfluenceFromFile(string fileName)
{
    unordered_map <ui, unordered_map <ui, double> > nodeNodeInfMap;

    ifstream uInfFile;
    uInfFile.open(fileName);

    string line;
    stringstream ss; 

    ui uNode, vNode, count;
    double infVal;

    if(uInfFile.is_open())
    {
        while(getline(uInfFile, line))
        {
            ss.clear();
            ss.str("");

            ss << line;

            ss >> count >> uNode;

            unordered_map <ui, double> tempMap;

            while(ss >> vNode >> infVal)
            {
                    tempMap[vNode] = infVal;
            }

            nodeNodeInfMap[uNode] = tempMap;
            tempMap.clear();
        }
    }
    else
    {
        cout << "Error opening file -- " << fileName << "\n";
    }
    return nodeNodeInfMap;
}
