# this code is for generating some random distributions sampled from Dirichlet...
# here we generate topic-word, topic-topic, and user-topic distributions 

from __future__ import division
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from collections import defaultdict
from collections import Counter

numTopics = 5
vocabsize = 500


def getUserBaseRates(fileName):

	tempUBR = {}

	f = open(fileName)

	for line in f:
		splts = line.strip().split()

		uid = int(splts[0].strip())
		ubr = float(splts[1].strip())

		tempUBR[uid] = ubr


	f.close()
	return tempUBR;


topicTopicVecFile = open("topicTopicDistribution.txt", "w")
hyperBeta = [0.01]*numTopics
topicTopicProbVectors = []
for x in range(numTopics):
	np.random.seed()
	
	localTopicTopicVec = list(np.random.dirichlet(hyperBeta))
	topicTopicProbVectors.append(localTopicTopicVec)

	writeStr = str(x) + ' ' + ' '.join([str(i) for i in localTopicTopicVec])

	topicTopicVecFile.write(writeStr + "\n")

topicTopicVecFile.close()


wordDistTopics = []
hyperAlpha = [0.1]*vocabsize
wordDistFile = open("topicWordDistribution.txt", "w")

for x in range(numTopics):
	np.random.seed()
	localWordDist = list(np.random.dirichlet(hyperAlpha))
	wordDistTopics.append(localWordDist)

	writeStr = str(x) + ' ' + ' '.join([str(i) for i in localWordDist]) 

	wordDistFile.write(writeStr + "\n")

wordDistFile.close()


print("Getting user base rates...")
userBaseRates = getUserBaseRates("sampleUsersBaseRate.txt")
print("Got user base rates", len(userBaseRates))


usersThatEmit = list(userBaseRates.keys())
print("users that emit -- ", len(usersThatEmit))


hyperGamma = [0.01]*numTopics
userTopicPrefVectors = defaultdict(list)
userTopicPrefFile = open("userTopicDistribution.txt", "w")

for x in usersThatEmit:
	np.random.seed()
	
	localUserTopicPref = list(np.random.dirichlet(hyperGamma))
	userTopicPrefVectors[x] = localUserTopicPref

	writeStr = str(x) + ' ' + ' '.join([str(i) for i in localUserTopicPref]) 

	userTopicPrefFile.write(writeStr + "\n")

userTopicPrefFile.close()