# HMHP model

The cpp file *hmhp_EstAll_GroupedWuv.cpp* is the code for the HMHP model. To run the model, the input files required are:

* events file
* documents file
* followers map 
* indegree and outdegree file for each node 
* candidate parents for each event and a file containing exponential time difference value for each candidate parent for each event (this can be  precomputed from the events file and followers map):

To compile the code:
```
g++ -std=c++11 -Wall -O2 -g -o hmhpModel hmhp_EstAll_GroupedWuv.cpp
```

To execute the code, **argv[1]** = *BURN-IN*, **argv[2]** = *Total No. Of Iterations*, **argv[3]** = *path to input files*, and **argv[4]** = *path to output files*. Code can be executed as follows:
```
./hmhpModel 200 301 inputFilesOurModel.txt outputFilesOurModel.txt
```

As we record the inferred parent and topic assignments only at every 10th iteration after the BURN-IN period, the total number of iterations must be greater than 10 and greater than BURN-IN period.

The code outputs various files as follows:

- parent assignment (recorded at every 10th iteration after BURN-IN)
- topic assignment (recorded at every 10th iteration after BURN-IN)
- avg parent assignment (probability of candidate parent event)
- grouped wuv values

After the execution of the above code and once all the files are in place, one can execute the script *evaluateHMHP.sh*

### The file format for the input files is as follows:

1. *Events File*: Each line of events file describes an event with 5 space separated values as follows:
```
event_time user_node parent_event_id topic_id level_info
```

2. *Documents File*: i-th line has a document corresponding to i-th event in the events file. Each document is set of space separated (int) word-ids.

3. *Followers Map*: Each line of the followers map is as follows (all the values are space separated):
```
NumOfFollowers userId followerId-1 followerId-2 ... followerId-M
```

4. *Indegree and Outdegree Files*: Each line of the both the files has two space separated values:
```
out/in-degree userId
```

5. *Candidate Parents File*: Each line contains the information of the candidate parent events for each event. Each line has space separated values as follows:
```
NumOfCandParents EventId CandParentId_1 CandParentId_2 ... CandParentId_100 
```

6. *Candidate Parents Exponential Time-difference file*: Each line contains the information of the exponential time difference between the candidate current event and the candidate parent event. Each line has space separated values as follows:
```
NumOfCandParents EventId ExpTimeDiffWithCandParentId_1 ExpTimeDiffWithCandParentId_2 ... ExpTimeDiffWithCandParentId_100
```

To get the candidate parent files, run the code *"getStoreTopKCandParents.cpp"*. The program takes the events file as input and outputs two files -- one for the candidate parents and other is the exponential time difference file. The program can be run as follows:
```
g++ -std=c++11 -Wall -O2 -g -o getTop100CandParents  getStoreTopKCandParents.cpp

./getTop100CandParents pathToEventsFile pathToCandidateParentFile pathToCandidateParentExpTimeDiffFile
```
