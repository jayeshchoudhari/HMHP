#include "namespace.h"
#include "dataIO.h"
#include "utilities.h"


using namespace std;

DataIO :: DataIO(string ipFilePathsFileName, string opFilePathsFileName)		// constructor
{
    inputFilePathsFileName = ipFilePathsFileName;
    outputFilePathsFileName = opFilePathsFileName;
    configInputFiles = getConfigInputOutputFileNames(inputFilePathsFileName);
    configOutputFiles = getConfigInputOutputFileNames(outputFilePathsFileName);
    populateDataInDataStructure(unordered_map<string, string> configInputFiles);
}


unordered_map<string, string> DataIO :: getConfigInputOutputFileNames(string fileName)
{
    ifstream inputOutputFiles;
	inputOutputFiles.open(fileName);

	unordered_map <string, string> localFilePathsObj;

	if(inputOutputFiles.is_open())
	{	
		string line, fileKey, fname;
		stringstream ss;

		while(getline(inputOutputFiles, line))
		{
			ss.clear();
			ss.str("");

			ss << line;
			ss >> fileKey >> fname;

			localFilePathsObj[fileKey] = fname;
		}
	}
	else
	{
		cout << "Issue opening the file containing the path to input files... Exiting here...\n";
		exit(0);
	}

	return localFilePathsObj;
}

int DataIO :: populateDataInDataStructure(unordered_map<string, string> configInputFiles)
{
    cout << "Reading Events...\n";
	allEvents = getSyntheticEventsFromFile(configInputFiles["eventsFile"]);
	cout << "Got all the events... " << allEvents.size() << "\n";

	allDocs = getSyntheticDocsFromFile(configInputFiles["docsFile"]);
	cout << "Got all the docs...\n";

    cout << "Creating hist of words for each doc...\n";
	createWordHistForAllDocs();
	cout << "Created hist of words for each doc...\n";

    cout << "Getting followers map\n"; 
	followersMap = readIntVectorMapFromFile(configInputFiles["mapFile"]);			// this would be required for the Wuv matrix.... 
	cout << "Got the followers map... " << followersMap.size() << "\n";

    cout << "Reading outdegree and indegree maps\n";
	outDegreeMap = getDegreeMap(configInputFiles["outDegreeFile"]);
	inDegreeMap = getDegreeMap(configInputFiles["inDegreeFile"]);
    cout << "Got the outdegree and indegree maps\n";

    // just to initialize alpha, beta, gamma, (Dirichlet prior) and sum of these vectors..
	readHyperParametersVectorsFromFile();	// for now not required...
	
    // populate all possible parent events...
	cout << "Getting possible parent events for all the events...\n";
    allPossibleParentEvents = populateParentEventsIds(configInputFiles["top100CandidateParentsFile"]);
    allPossibleParentEventsExponentials = populateParentEventsExpKernelValue(configInputFiles["top100CandidateParentsExpFile"]);
	cout << "Got candidate parent events and their exp latency in time..." << allPossibleParentEvents.size() << " " << allPossibleParentEventsExponentials.size() << "\n";

    return 0;
}

// can be moved to utilities if documents are passed...
int DataIO :: createWordHistForAllDocs()
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

// can be moved to utilities... 
unordered_map <ui, ui> DataIO :: getHistOfWordsOverWordsFromDoc(vector<ui> doc)
{
	unordered_map <ui, ui> wordCounts;

	for(unsigned int i = 0; i < doc.size(); i++)
	{
		li word = doc[i];

		wordCounts[word] += 1;
	}

	return wordCounts;
}


vector < vector <ui> > DataIO :: getSyntheticDocsFromFile(string fileName)
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


vector < vector <ui> > DataIO :: readIntVectorMapFromFile(string fileName)
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

unordered_map<int, int> DataIO :: getDegreeMap(string degreeMapFile)
{
	unordered_map<int, int> degreeMap;

	ifstream dMFile;
	dMFile.open(degreeMapFile);

	string line;
	stringstream ss;

	int deg, uNode;

	if(dMFile.is_open())
	{
		while(getline(dMFile, line))
		{
			ss.clear();
			ss.str("");

			ss << line;

			ss >> deg >> uNode;

			degreeMap[uNode] = deg;
		}
	}

	return degreeMap;
}

int DataIO :: readHyperParametersVectorsFromFile()
{
	// hyperBeta = readDoubleVectorFromFile("hyperBetaValue.txt");
	// hyperGamma = readDoubleVectorFromFile("hyperGammaValue.txt");
	// hyperAlpha = readDoubleVectorFromFile("hyperAlphaValue.txt");


	// topicTopicProbVector = readMultipleDoubleVectorsFromFile("topicTopicProbVectors.txt");
	// userTopicPrefVector = readMultipleDoubleVectorsFromFile("nodeTopicProbVectors.txt");
	// topicWordProbVector = readMultipleDoubleVectorsFromFile("topicWordProbVectors.txt");

	vector<double> tempBeta(numTopics, 0.1);
	hyperBeta = tempBeta;

	vector<double> tempGamma(numTopics, 0.01);
	hyperGamma = tempGamma;

	vector<double> tempAlpha(vocabSize, 0.1);
	hyperAlpha = tempAlpha;

	sumAlpha = accumulate(hyperAlpha.begin(), hyperAlpha.end(), 0.0);
	sumBeta = accumulate(hyperBeta.begin(), hyperBeta.end(), 0.0);
	sumGamma = accumulate(hyperGamma.begin(), hyperGamma.end(), 0.0);
	return 0;
}

vector <vector <ui> > DataIO :: populateParentEventsIds(string parentsFile)
{
    vector <vector <ui> > localAllPossibleParentEvents;
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

			localAllPossibleParentEvents.push_back(parEvents);
			parEvents.clear();
		}

		possParFile.close();
	}
	else
	{
		cout << "Cannot open the possible parent events file\n";
		exit(0);
	}

    return localAllPossibleParentEvents;
}

vector <vector <double> > DataIO :: populateParentEventsExpKernelValue(string parentExpFile)
{
    vector <vector <double> > localAllPossibleParentEventsExponentials;
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
	
	return localAllPossibleParentEventsExponentials;
}