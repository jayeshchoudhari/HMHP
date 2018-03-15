from __future__ import division
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from collections import defaultdict
from collections import Counter

lastLevel = 3

eventsFileName = "../events_semisyn_sample_gamma_1M.txt"
followersFileName = "../../mapped_users_followers_restricted_tweets_top_5k_hashes.txt"
outDegreeFileName = "../outdegree_nodeid.txt"
inDegreeFileName = "../indegree_nodeid.txt"
# wuvOutFile = "wuvEstFile_1Lakh.txt"
wuvOutFile = "wuvEstFile_fromLastIterationParent.txt"


logBase = 1.5
baseAlpha = 0.01
baseBeta = 1
nonZeroNuSumElements = defaultdict(list)
staticGid = "0_0"

#  get assigned parent for each event....
def getLevelInfo():

	topParentsFile = open("./parentAssignmentLastIteration.txt")

	levelOCounts = 0
	# levelVal = []
	localParentInfo = []

	# not required for now...
	becameParent = [-1]*1000001

	maxLevel = 0
	linenum = 0

	for line in topParentsFile:
		flds = line.strip().split()
		eid = int(flds[0].strip())
		pid = int(flds[1].strip())

		if pid > -1:
			nearParentEvent = int(flds[1].strip())
			# parLevel = levelVal[nearParentEvent]
			# levelVal.append(parLevel + 1)
			localParentInfo.append(nearParentEvent)
		else:
			levelOCounts += 1
			# levelVal.append(0)
			localParentInfo.append(-1)

		linenum += 1

	
	return [becameParent, localParentInfo]


# read the original set of events
# also get the count of events for each node from the first 15% of the events...
def getAllEvents(fileName):

	localAllEvents = []
	localNodeEventsCountMap = defaultdict(int)

	eventsFile = open(fileName)

	totalSum = 0
	actualSum = 0
	linenum = 0
	contrib = 0
	for line in eventsFile:
		singleEvent = []
		flds = line.strip().split()
		localAllEvents.append([float(flds[0].strip()), int(flds[1].strip()), int(flds[2].strip()), int(flds[3].strip())])

		if linenum <= 150000:
			# increase the count of event for a node...
			localNodeEventsCountMap[int(flds[1].strip())] += 1
			totalSum += 1
			# contrib += 1
		
		linenum += 1
	
	eventsFile.close()
	
	# print "contrib = ", contrib
	print "total sum = ",  totalSum
	print "actual sum = ", actualSum
	return [localAllEvents, localNodeEventsCountMap]

'''
def getNodeEventsCount():

	localNodeEventsCountMap = defaultdict(int)

	for i in range(len(allEvents)):

		singleEvent = allEvents[i]
		uNode = singleEvent[1]
		localNodeEventsCountMap[uNode] += 1

	return localNodeEventsCountMap
'''

# read the outdegree and indegree file
# outdegree and indegree file format -- each line has the two values:
# outdegree/indegree nodeid..
def getOutInDegreeValues(fileName):

	degreeFile = open(fileName)

	nodeDegrees = defaultdict(int)

	for line in degreeFile:
		flds = line.strip().split()
		outVal = int(flds[0].strip())
		nodeId = int(flds[1].strip())

		nodeDegrees[nodeId] = outVal

	degreeFile.close()

	return nodeDegrees

# use the map
# using the map get the Nu counts...
# go to each edge and increase the Nu count for the group by the number of events (in the first 15%) of the parent node...
def getGroupSourceDenomSum(fileName):

	global nonZeroNuSumElements

	localGroupSourceDenomSum = defaultdict(int)

	follFile = open(fileName)

	for line in follFile:
		flds = line.strip().split()
		count = int(flds[0].strip())

		uNode = int(flds[1].strip())
		outDegVal = outDegreeMap[uNode]
		outDegGid = int(math.floor(math.log(outDegVal, logBase)))	
		outGid = str(outDegGid)

		for j in range(2, len(flds)):

			vNode = int(flds[j].strip())
			inDegVal = inDegreeMap[vNode]
			inDegGid = int(math.floor(math.log(inDegVal, logBase)))	
			inGid = str(inDegGid)

			gid = outGid + "_" +  inGid

			# if uNode is not present in the nodeEventsCount dict then the count is zero...
			# localGroupSourceDenomSum[gid] += nodeEventsCount[uNode]
			localGroupSourceDenomSum[staticGid] += nodeEventsCount[uNode]

			# just to get some stats...
			if nodeEventsCount[uNode] > 0:
				nonZeroNuSumElements[staticGid].append(nodeEventsCount[uNode])

	follFile.close()

	return localGroupSourceDenomSum

#  using the assigned parents get the Nuv counts (for each edge)...
def getNodeNodeCountMap():

	localNodeNodeCount = defaultdict(lambda: defaultdict(int))

	for i in range(len(allEvents)):

		singleEvent = allEvents[i]

		childNode = singleEvent[1]
		# pid = singleEvent[2]
		pid = parentInfo[i]
		# parentInfo contains the information of the assigned parents 

		if pid > -1:
			# parNode = allEvents[pid][1]
			# get the node of the parent event...
			parNode = allEvents[pid][1]
			localNodeNodeCount[parNode][childNode] += 1

	return localNodeNodeCount

# get the group based Nuv value...
def getGroupSumNumePerGid():

	localGroupNumeSum = defaultdict(int)

	for uNode in nodeNodeCountMap:

		outDegVal = outDegreeMap[uNode]
		outDegGid = int(math.floor(math.log(outDegVal, logBase)))	
		outGid = str(outDegGid)

		for vNode in nodeNodeCountMap[uNode]:

			inDegVal = inDegreeMap[vNode]
			inDegGid = int(math.floor(math.log(inDegVal, logBase)))	
			inGid = str(inDegGid)

			gid = outGid + "_" +  inGid

			# localGroupNumeSum[gid] += nodeNodeCountMap[uNode][vNode]
			localGroupNumeSum[staticGid] += nodeNodeCountMap[uNode][vNode]
	
	return localGroupNumeSum


# estimate Wuv using the group based Nuv and Nu values...
def estimateWuv():

	outFile = open(wuvOutFile, "w")

	for gid in groupNumeSum:

		# wuvVal = (groupNumeSum[gid] + baseAlpha) / (groupSourceDenomSum[gid] + baseBeta)
		wuvVal = (groupNumeSum[staticGid] + baseAlpha) / (groupSourceDenomSum[staticGid] + baseBeta)

		# outFile.write(gid + " " + str(wuvVal) + " " + str(groupNumeSum[gid]) + " " + str(groupSourceDenomSum[gid]) + "\n")
		outFile.write(staticGid + " " + str(wuvVal) + " " + str(groupNumeSum[staticGid]) + " " + str(groupSourceDenomSum[staticGid]) + "\n")

	outFile.close()


# get the assigned parents info...
eventsThatBecameParent, parentInfo = getLevelInfo()
print "Got parent information...", len(parentInfo)

# get the original events and the count of the events for each node from first 15% of the events....
allEvents, nodeEventsCount = getAllEvents(eventsFileName)
print("Got all the events -- ", len(allEvents))

# get outdegree and indegree for each node...
outDegreeMap = getOutInDegreeValues(outDegreeFileName)
inDegreeMap = getOutInDegreeValues(inDegreeFileName)
print("Got In and out degrees -- ", len(outDegreeMap), len(inDegreeMap))

# get the Nu values for each group
print("Reading followers map...")
groupSourceDenomSum = getGroupSourceDenomSum(followersFileName)
print("Got Group denom sum")

#  get the Nuv for each edge 
print("getting node node counts")
nodeNodeCountMap = getNodeNodeCountMap()
# print(nodeNodeCountMap)

# get group based Nuv
groupNumeSum = getGroupSumNumePerGid()

# estimate Wuv...
print("estimating wuv")
estimateWuv()
