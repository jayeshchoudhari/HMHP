#include "./include/namespace.h"
#include "./include/dataIO.h"
#include "./include/gibbsSampler.h"
#include "./include/utilities.h"
#include <bits/stdc++.h>

using namespace std;

int main(int argc, char *argv[])
{
	int numOfTopics = 5;		// num of Topics
	ui maxNumOfNodes = 1000002;	// num of nodes in user graph
	int BURN_IN = 200;			// default values
	ui ITERATIONS = 301;		// default values
	
	string inputFilePaths, outputFilePaths;

	if(argc == 5)
	{
		BURN_IN = atoi(argv[1]);
		ITERATIONS = atoi(argv[2]);

		inputFilePaths = argv[3];
		outputFilePaths = argv[4];
        
		cout << "Will be running sampler for BURN_IN = " << BURN_IN << " And ITERATIONS = " << ITERATIONS << "\n";
	}
	else
	{
		cout << "Mismatch in the command line arguments... \n argv[1] = BURN_IN, argv[2] = TotalIterations(including BURN_IN), argv[3] = PathToInputFiles, argv[4] = PathToOutputFiles\n";
		exit(0);
	}
    
	DataIO dataIOObj(inputFilePaths, outputFilePaths, numOfTopics, maxNumOfNodes);
	cout << "Read the data...\n";

	GibbsSampler GBS(dataIOObj);
	cout << "Done with Sampler Initializing...\n";

	GBS.IterativeSampler(dataIOObj, ITERATIONS, BURN_IN);

	ofstream llFile;
	// llFile.open("estAllLogLikelihood.txt");
	llFile.open(dataIOObj.configOutputFiles["logLikelihoodFile"]);



    return 0;
}