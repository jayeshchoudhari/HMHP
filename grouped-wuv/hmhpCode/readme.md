
The cpp file hmhp_EstAll_GroupedWuv.cpp is the code for the HMHP model. To run the model, the input files required are events, documents, followers map, indegree and outdegree for each node, candidate parents for each event and exponential time difference value for each candidate parent for each event (this could be precomputed from the events file and followers map.):

To compile the code:
	g++ -std=c++11 -Wall -O2 -g -o hmhpModel hmhp_EstAll_GroupedWuv.cpp

To execute the code, argv[1] = BURN-IN, argv[2] = Total No. Of Iterations, argv[3] = path to input files, and argv[4] = path to output files. Code can be executed as follows:

	./hmhpModel 200 301 inputFilesOurModel.txt outputFilesOurModel.txt


As we record the inferred parent and topic assignments only at every 10th iteration after the BURN-IN period, the total number of iterations must be greater than 10 and greater than BURN-IN period.

The code outputs various files as follows:
	-- parent assignment (recorded at every 10th iteration after BURN-IN)
	-- topic assignment (recorded at every 10th iteration after BURN-IN)
	-- avg parent assignment (probability of candidate parent event)
	-- grouped wuv values

After the execution of the above code and once all the files are in place, one can execute the script "evaluateHMHP.sh"
