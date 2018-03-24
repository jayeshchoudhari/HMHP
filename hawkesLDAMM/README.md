# Hawkes + LDA Model

The model implemented here performs *topic modeling*, *parent identification*, and *network reconstruction* (each of the) tasks independently.

*Topic Modeling* task is performed using the LDA-Mixture Model. Here each document has a single topic.

*Parent Identification* task is done using the Hawkes process. It turns out that *top100ParentsFile.txt* files in the *inputFiles* folder is exactly the output of the *parent identification* task. 


*Network Reconstruction*: Once the *parent identification* task is done, then the network reconstruction task is performed using *estimateWuvFromNearestParent_withoutPrior.py*. The edges of the network are grouped on the basis of the outdegree of source node and indegree of the destination node for the *network reconstruction* task.