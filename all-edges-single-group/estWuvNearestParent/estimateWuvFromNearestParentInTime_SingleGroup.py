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
wuvOutFile = "wuvEstFile_inducedLevel.txt"


logBase = 1.5
baseAlpha = 0.01
baseBeta = 1
nonZeroNuSumElements = defaultdict(list)
staticGid = "0_0"

def getLevelInfo():
	topParentsFile = open("../top100CandParents_sample_events_1L_time_15.txt")

	levelOCounts = 0
	# levelVal = []
	becameParent = [-1]*100001
	localParentInfo = []

	maxLevel = 0
	linenum = 0

	for line in topParentsFile:
		flds = line.strip().split()
		count = int(flds[0].strip())
		eid = int(flds[1].strip())

		if count > 0:
			nearParentEvent = int(flds[2].strip())
			# parLevel = levelVal[nearParentEvent]
			# levelVal.append(parLevel + 1)
			# becameParent[nearParentEvent] = 1
			localParentInfo.append(nearParentEvent)
		else:
			levelOCounts += 1
			# levelVal.append(0)
			localParentInfo.append(-1)

		linenum += 1

	
	# leafEvents = list(filter(lambda x: x == -1, becameParent))
	# print "# leaf Events = ", len(leafEvents)

	return [becameParent, localParentInfo]


# read the events
def getAllEvents(fileName):

	localAllEvents = []
	localNodeEventsCountMap = defaultdict(int)

	eventsFile = open(fileName)

	# leafEvents = list(filter(lambda x: x == -1, eventsThatBecameParent))
	# print "# leaf Events = ", len(leafEvents)

	totalSum = 0
	actualSum = 0
	linenum = 0
	contrib = 0
	for line in eventsFile:
		singleEvent = []
		flds = line.strip().split()
		localAllEvents.append([float(flds[0].strip()), int(flds[1].strip()), int(flds[2].strip()), int(flds[3].strip())])

		# localNodeEventsCountMap[int(flds[1].strip())] += 1		
		# if int(flds[4].strip()) != lastLevel:
			# actualSum += 1

		# if eventsThatBecameParent[linenum] != -1:
        # if linenum <= 25000 and eventsThatBecameParent[linenum] != -1:
		if linenum <= 150000:
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

			# localGroupSourceDenomSum[gid] += nodeEventsCount[uNode]
			localGroupSourceDenomSum[staticGid] += nodeEventsCount[uNode]

			if nodeEventsCount[uNode] > 0:
				nonZeroNuSumElements[gid].append(nodeEventsCount[uNode])

	follFile.close()

	return localGroupSourceDenomSum


def getNodeNodeCountMap():

	localNodeNodeCount = defaultdict(lambda: defaultdict(int))

	for i in range(len(allEvents)):

		singleEvent = allEvents[i]

		childNode = singleEvent[1]
		# pid = singleEvent[2]
		pid = parentInfo[i]

		if pid > -1:
			# parNode = allEvents[pid][1]
			parNode = allEvents[pid][1]
			localNodeNodeCount[parNode][childNode] += 1

	return localNodeNodeCount


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


def estimateWuv():

	outFile = open(wuvOutFile, "w")

	for gid in groupNumeSum:

		# wuvVal = (groupNumeSum[gid] + baseAlpha) / (groupSourceDenomSum[gid] + baseBeta)
		wuvVal = (groupNumeSum[staticGid] + baseAlpha) / (groupSourceDenomSum[staticGid] + baseBeta)

		outFile.write(staticGid + " " + str(wuvVal) + " " + str(groupNumeSum[staticGid]) + " " + str(groupSourceDenomSum[staticGid]) + "\n")

	outFile.close()



eventsThatBecameParent, parentInfo = getLevelInfo()
print "Got parent information...", len(parentInfo)

allEvents, nodeEventsCount = getAllEvents(eventsFileName)
print("Got all the events -- ", len(allEvents))

outDegreeMap = getOutInDegreeValues(outDegreeFileName)
inDegreeMap = getOutInDegreeValues(inDegreeFileName)
print("Got In and out degrees -- ", len(outDegreeMap), len(inDegreeMap))


print("Reading followers map...")
groupSourceDenomSum = getGroupSourceDenomSum(followersFileName)
print("Got Group denom sum")

print("getting node node counts")
nodeNodeCountMap = getNodeNodeCountMap()
# print(nodeNodeCountMap)

groupNumeSum = getGroupSumNumePerGid()

print("estimating wuv")
estimateWuv()
