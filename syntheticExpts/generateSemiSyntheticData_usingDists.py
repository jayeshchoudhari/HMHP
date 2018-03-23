from __future__ import division
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from collections import defaultdict
from collections import Counter

#################### Global Part ###########################

Tmax = int(sys.argv[1])
alphaMultiplier = float(sys.argv[2])
seedVal = int(sys.argv[3])
np.random.seed(seedVal)

# Tmax = 100
# alphaMultiplier = 1
# seedVal = 1234
# np.random.seed(seedVal)
# eventsToBeGenerated = 1000000
eventsToBeGenerated = int(sys.argv[4])

numTopics = int(sys.argv[5])
# numTopics = 100

vocabsize = int(sys.argv[6])
# vocabsize = 33900

maxLevelsToBeGenerated = 25

# Output Files...
eventsFileName = "test_level_events_semisyn_sample.txt"
docsFileName = "test_level_docs_semisyn_sample.txt"

allSyntEventsFileName = eventsFileName;
allSyntDocsFileName = docsFileName

# Input Files...
# userUserInfluenceFileName = "../semi_synthetic_expt_diag_grouped_wuv_all_edges_ssegwall/semiSyntheticWuvValues_top_5k_hashes_Grouped_All.txt"
# userBaseRatesFileName = "calculatedUserBaseRate_FullDay.txt"
# userTopicPrefVectorsFileName = "userTopicDistribution.txt"
# topicTopicProbVectorsFileName = "topicTopicDistribution_1M.txt"
# wordDistTopicsFileName = "topicWordDistribution_1M.txt"

userUserInfluenceFileName = "./sampleUserUserInfluence.txt"
userBaseRatesFileName = "./sampleUsersBaseRate.txt"
userTopicPrefVectorsFileName = "./userTopicDistribution.txt"
topicTopicProbVectorsFileName = "./topicTopicDistribution.txt"
wordDistTopicsFileName = "./topicWordDistribution.txt"


# numTopics = 100
topics = [i for i in range(0,numTopics)]

docsize = 10

rayleighSigma = 3

alphaArrivsWuv = defaultdict(list)
alphaArrivsInterArrivs = defaultdict(list)
meanWuvs = defaultdict(list)

#################### Global Part ###########################

#################### Functions ########################

def get_lambda_for_level_l(startTime, x, W_uv):
	
	timeDiff = x - startTime
	# func_value = W_uv * math.exp(-timeDiff)
	func_value = W_uv * math.exp(- alphaMultiplier * timeDiff)

	return func_value


def get_lambda_for_level_l_rayleigh(startTime, x, W_uv):

	timeDiff = x - startTime
	# func_value = W_uv * math.exp(-timeDiff)
	xBySigmaSqr = timeDiff*1.0/(rayleighSigma**2)

	xSqrBy2SigmaSqr = (timeDiff**2)*1.0 / (2*(rayleighSigma**2))

	func_value = W_uv * xBySigmaSqr * math.exp(-xSqrBy2SigmaSqr)

	return func_value	


def get_nhpp_timestamps_for_level_l(Tmax, lambdaUpperBound, t_e, W_uv):

	T, t0, n, m = Tmax, t_e, 0, 0
	lambda_constant = lambdaUpperBound 
	sm = t_e
	# homoTime = []
	inhomoTime = []
	# allDs = []

	while sm < T:
		u = np.random.uniform(0,1)
		w = -(np.log(u) / lambda_constant)						# so that w ~ exp(lambda_constant)
		sm = sm + w 
		probVal = (get_lambda_for_level_l(t_e, sm, W_uv))/lambda_constant
        # probVal = (get_lambda_for_level_l_rayleigh(t_e, sm, W_uv))/lambda_constant
		if sm < T and probVal >= 1e-7:
			# homoTime.append(sm)									# sm are the points in the homo. PP
			d = np.random.uniform(0,1)						
			# print(sm, d, math.exp(-(sm - t_e)))
			if d <= probVal:
				# allDs.append(d)
				inhomoTime.append(sm)							# inhomoTime are the points in inhomo. PP
				# n += 1
			# m += 1

		else:
			break

	# print (inhomoTime)
	# plot_for_confirmation(homoTime, inhomoTime, allDs)
	return inhomoTime


def get_homogenous_pp_timestamps(Tmax, lambda_rate):

	T, t0 = Tmax, 0
	hppTimeStamps = []
	prev_ti = t0

	while prev_ti <= Tmax:
		u = np.random.uniform(0,1)
		x = prev_ti - (np.log(u)/lambda_rate)
		if x <= Tmax:
			hppTimeStamps.append(x)
		prev_ti = x

	return hppTimeStamps


# def generate_synthetic_events(nodes, Tmax, maxConstLambda, lamdaFunc, W, muVal):
def generate_synthetic_events():

	# numNodes = len(nodes)
	noparent = '-1'
	allEvents = []
	# will be a sorted list of events
	# where each event is a list [time, node, parent, topic]
	# topic to each doc/event will be appended later...
	eventsWithId = {}

	level0Events = []

	# Generate Level-0 events

	# get a random set of 100 nodes... and generate L0-events only for those nodes...
	# randomNodes = np.random.choice(usersThatEmit, 100)
	# for node in range(numNodes):
	# for node in nodes:
    # for node in usersThatEmit:
    # for node in randomNodes:
	for node in usersThatEmit:

		# utils.get_nhpp_timestamps(Tmax, lambdaFunction, lambdaUpperBound) 
		# where lambdaUpperBound is a constant lambda upper bounding the rate or intensity function
		# currently we are assuming that is known for each node
		# nodeL0Timestamps = utils.get_nhpp_timestamps(Tmax, lamdaFunc, maxConstLambda)
		userMuVal = userBaseRates[node]
		# userMuVal = 0.02
		nodeL0Timestamps = get_homogenous_pp_timestamps(Tmax, userMuVal)

		# each event is a list [tn, cn, pn, eta_n] = [time, node, parent, topic]...
		# topic for each event will be assigned later...
		# we use -1 as parent for level0 events as first event has eventId 0
 		if userTopicPrefVectors.get(node, -1) != -1:

			for ts in nodeL0Timestamps:
				probVec = np.array(userTopicPrefVectors[node])
				probVec /= probVec.sum().astype(float)
				userTopic = np.random.choice(topics, 1, True, probVec)[0]
		 		level0Events.append([ts, node, noparent, userTopic, 0])


	# sort the events with timestamps
	level0Events.sort(key = lambda row:row[0])

	print("Len level 0 events = ", len(level0Events))
	for event in level0Events:
		allEvents.append(event)

	# print(len(allEvents))
	
	eventsCount = len(allEvents)
	# Generate level-l events
	prevLevelEvents = level0Events
	currentLevelEvents = []

	# print("Done with level 0 events...")
	level = 1
	nomoreEventsFlag = 1

	while len(prevLevelEvents) > 0:
		# print("No. Prev level Events -- ", len(prevLevelEvents))
		currentLevelEvents[:] = []
		for event in prevLevelEvents:
			u_ts = event[0]
			u = event[1]
			parentEventId = str(u) + "_" + str(u_ts)
			parentTopic = event[3]

			upresent = W.get(u, -1)
			# if node u is present in the followers map or if u has followers then go for generating events from the neighbors of u...
			# if type(upresent) == dict:
			if upresent != -1:

				for node in W[u]:

					W_uv = W[u][node]
					# nhppTimestampsLevelL = utils.get_nhpp_timestamps_for_level_l(Tmax, 1, u_ts, W_uv)
					# nhppTimestampsLevelL = utils.get_nhpp_timestamps_for_level_l(Tmax, lambda_constant, u_ts, W_uv)

					if W_uv > 0:

						nhppTimestampsLevelL = get_nhpp_timestamps_for_level_l(Tmax, W_uv, u_ts, W_uv)						# lambda_constant = W_uv
						# tempEvents = []

						srcDest = str(u) + "_" + str(node)
						inteArrivSum = 0
						for ts in nhppTimestampsLevelL:
							probVec = np.array(topicTopicProbVectors[parentTopic])

							if probVec.sum() == 0:
								eta_n = np.random.choice(topics, 1)[0]
							else:
								probVec /= probVec.sum().astype(float)
								eta_n = np.random.choice(topics, 1, True, probVec)[0]
							
							currentLevelEvents.append([ts, node, parentEventId, eta_n, level])
							eventsCount += 1
						# inteArrivSum += (ts - u_ts)

					'''
					# meanWuvs[srcDest].append(len(nhppTimestampsLevelL))
					if inteArrivSum > 0:
						alphaArrivsWuv[srcDest].append(W_uv)
						alphaArrivsInterArrivs[srcDest].append(inteArrivSum)
						# alphaArrivs[srcDest].append(math.sqrt(W_uv/Tmax))
					'''
					# currentLevelEvents = currentLevelEvents + tempEvents			# concatenation of two lists...

                    # if eventsCount % 100000 == 0:
                        # print(eventsCount)

					if eventsCount > eventsToBeGenerated:
						nomoreEventsFlag = 0
						break


			if nomoreEventsFlag == 0:
				break

		# append new events to the list of existing global (aggregate) events...
		# could just add two lists...
		if len(currentLevelEvents) > 0:
			allEvents = allEvents + currentLevelEvents
			prevLevelEvents = list(currentLevelEvents)

		else:
			prevLevelEvents[:] = []

		print(len(allEvents))
		print("Done with level", level, "events") 

		level += 1

		if level > maxLevelsToBeGenerated:
			print("done with required levels... breaking here....")
			break

		if len(allEvents) > eventsToBeGenerated:
			print("Events more than 1000000.. Breaking here")
			break

	allEvents.sort(key = lambda row:row[0])

	# print("Got all the events... lets have the parent structure in place...")

	# this will map the uniqEventId with its index in the sorted list of events....
	uniqEventIdSortedIndex = {}
	for eid in range(len(allEvents)):
		ts = allEvents[eid][0]
		node = allEvents[eid][1]
		uniqEventId = str(node) + "_" + str(ts)
		uniqEventIdSortedIndex[uniqEventId] = eid

	# print("Adding the required parent Index to the events...")
	# now for each event, if its parent is some uniqEventId, we have change that with the event index....
	for eid in range(len(allEvents)):
		uniqParentId = allEvents[eid][2]
		if uniqParentId != '-1':
			indexOfParent = uniqEventIdSortedIndex[uniqParentId]
			allEvents[eid][2] = indexOfParent
		else:
			allEvents[eid][2] = -1

	print("numEvents = ", len(allEvents))

	return allEvents


def writeOnlyEventsToFile():

	print("Writing events to file....\n")

	allEventsFile = open(allSyntEventsFileName, 'w')
	# Sample topic for each event based on the parent
	for event in range(0, len(allSyntheticEvents)): 
        # allSyntheticEvents[event].append(-1)

		# write synthetic data to file...
		eventStr = ' '.join([str(ele) for ele in allSyntheticEvents[event]])
		eventStr = eventStr + "\n"
		allEventsFile.write(eventStr)

	allEventsFile.close()


# def get_assigned_topics_events(nodes, topics, allSyntheticEvents, numTopics, topic_topic, user_topic):
def get_assigned_topics_events():

	print("Assigning Topics to events...");
	# numNodes = len(nodes)
	
	allEventsFile = open(allSyntEventsFileName, 'w')

	# Sample topic for each event based on the parent
	for event in range(0, len(allSyntheticEvents)): 
		tn, un, pn = allSyntheticEvents[event][0], allSyntheticEvents[event][1], allSyntheticEvents[event][2]
		# pn here gives the eventId (index) of the parent event if the parent event is not -1
		if pn == -1:
		# if there is no parent.... sample from the user topic preferences...
			try:
				eta_n = np.random.choice(topics, 1, True, userTopicPrefVectors[un])[0]
			except:
				print("Key Error", un)
			#returns an array, so take the first index
		else:
		# if a parent exists...
		# pn is parent eventId (index)...
			# parentNode = allSyntheticEvents[pn][1]
			parentTopic = allSyntheticEvents[pn][3]
			eta_n = np.random.choice(topics, 1, True, topicTopicProbVectors[parentTopic])[0]
			#returns an array, so take the first index
			
		allSyntheticEvents[event].append(eta_n)

		# write synthetic data to file...
		eventStr = ' '.join([str(ele) for ele in allSyntheticEvents[event]])
		eventStr = eventStr + "\n"
		allEventsFile.write(eventStr)

	allEventsFile.close()

	# print(allSyntheticEvents)
	return allSyntheticEvents


# def generate_synthetic_docs(allSyntDocsFileName, allSyntheticEventsTopicsAssigned, vocabsize, numTopics, docsize):
def generate_synthetic_docs(allSyntDocsFileName):

	print("Generating Documents...")
	docsize = np.random.poisson(7)

	vocabulary = [i for i in range(0,vocabsize)]
	# Generate word distribution for each topic
	# hyperAlpha = [np.random.uniform(0,1) for k in range(0,vocabsize)]

	# Generate words for all the docs and write to file
	docWords = []

	allDocsFile = open(allSyntDocsFileName, 'w')

    # for eachdoc in allSyntheticEventsTopicsAssigned:
	for eachdoc in allSyntheticEvents:
		docsize = np.random.poisson(7)
		wordsInDoc = np.random.choice(vocabulary, docsize, True, wordDistTopics[eachdoc[3]])
		# draw doc size number of words
		docWords.append(wordsInDoc)
		docStr = ' '.join([str(word) for word in wordsInDoc])
		docStr = docStr + "\n"
		allDocsFile.write(docStr)
        # append event id

	allDocsFile.close()
	# print(docWords)
	return docWords



def getNodeList(fileName):

	nodesFile = open(fileName)

	tempNodeList = []

	for line in nodesFile:
		tempNodeList.append(int(line.strip()))

	return tempNodeList;


def getUserUserInfluence(fileName):

	wuvfile = open(fileName)

	tempWuv = {}

	for line in wuvfile:
		splitLine = line.split()
		count = int(splitLine[0].strip())
		uNode = int(splitLine[1].strip())

		tempWuv[uNode] = {}
		# followers[uNode] = []

		j = 2
		while j < len(splitLine):
			vNode = int(splitLine[j].strip())
			uvInf = float(splitLine[j+1].strip())
						
			tempWuv[uNode][vNode] = uvInf
			j = j + 2


	wuvfile.close()
	return tempWuv


def getFollowers(fileName):

	tempFollowers = {}

	f = open(fileName)

	i = 0
	for line in f:
		splts = line.strip().split()
		count = int(splts[0].strip())
		uid = int(splts[1].strip())

		uidFoll = []

		j = 2
		while j < len(splts):
			uidFoll.append(int(splts[j].strip()))
			j += 1


		tempFollowers[uid] = list(uidFoll)

		# print(i)
		# i += 1

	f.close()
	return tempFollowers


def getUserTopicPrefVectors(fileName):

	tempUserTopicPrefVector = {}

	f = open(fileName)

	for line in f:
		splts = line.strip().split()

		uid = int(splts[0].strip())

		tempVec = []
		sumVec = 0
		for j in range(1,len(splts)):
			ele = float(splts[j].strip())
			# sumVec += ele
			tempVec.append(ele)

		tempVec = np.array(tempVec)
		tempVec /= tempVec.sum()

		tempUserTopicPrefVector[uid] = list(tempVec)

	f.close()
	return tempUserTopicPrefVector


def getTopicWordProbVectors(fileName):

	tempWordDistVectors = []

	f = open(fileName)

	for line in f:
		splts = line.strip().split()

		tid = int(splts[0].strip())

		tempVec = []
		for j in range(1,len(splts)):
			tempVec.append(float(splts[j].strip()))

		tempVec = np.array(tempVec)
		tempVec /= tempVec.sum()

		tempWordDistVectors.append(list(tempVec))


	f.close()
	return tempWordDistVectors


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


def getTopicTopicProbVectors(fileName):

	tempTopicTopicVec = []

	f = open(fileName)

	for line in f:
		splts = line.strip().split()

		topId = int(splts[0])
		topicVec = []

		j = 1
		while j < len(splts):
			topicVec.append(float(splts[j]))
			j = j + 1

		topicVec = np.array(topicVec)
		topicVec /= topicVec.sum()

		tempTopicTopicVec.append(list(topicVec))

	return tempTopicTopicVec;




#################### Functions ###########################


############################## 
####### global part ##########
############################## 

# nodes = getNodeList(nodeListFileName)
# print("Got the list of nodes...", len(nodes))


'''
print("Getting followers map...")
followers = getFollowers(followersFileName)
print("Got the followers map...", len(followers)) 
'''

print("Getting user user influence map...")
W = getUserUserInfluence(userUserInfluenceFileName)
print("Got user-user influence matrix", len(W))

nodes = list(W.keys())

# numNodes = len(nodes)

print("Getting user base rates...")
userBaseRates = getUserBaseRates(userBaseRatesFileName)
print("Got user base rates", len(userBaseRates))


usersThatEmit = list(userBaseRates.keys())
print("users that emit -- ", len(usersThatEmit))


print("Getting Topic Topic Vectors...")
topicTopicProbVectors = getTopicTopicProbVectors(topicTopicProbVectorsFileName)

'''
topicTopicVecFile = open("topicTopicInteraction.txt", "w")
hyperBeta = [0.01]*numTopics
topicTopicProbVectors = []
for x in range(numTopics):
	np.random.seed()
	
	localTopicTopicVec = list(np.random.dirichlet(hyperBeta))
	topicTopicProbVectors.append(localTopicTopicVec)

	writeStr = str(x) + ' ' + ' '.join([str(i) for i in localTopicTopicVec])

	topicTopicVecFile.write(writeStr + "\n")

topicTopicVecFile.close()
'''
print("Got topic topic prob vectors", len(topicTopicProbVectors))


####################################

print("Getting word dist for each topic")
wordDistTopics = getTopicWordProbVectors(wordDistTopicsFileName)


'''
wordDistTopics = []
hyperAlpha = [0.1]*vocabsize
wordDistFile = open("OrigWordDist.txt", "w")

for x in range(numTopics):
	np.random.seed()
	localWordDist = list(np.random.dirichlet(hyperAlpha))
	wordDistTopics.append(localWordDist)

	writeStr = str(x) + ' ' + ' '.join([str(i) for i in localWordDist]) 

	wordDistFile.write(writeStr + "\n")

wordDistFile.close()
'''
print("Got word dist for each topic", len(wordDistTopics))

##########################


print("Getting topic dist for each user")
userTopicPrefVectors = getUserTopicPrefVectors(userTopicPrefVectorsFileName)
'''
hyperGamma = [0.01]*numTopics
userTopicPrefVectors = defaultdict(list)
userTopicPrefFile = open("userTopicPreference_level_restricted.txt", "w")

for x in usersThatEmit:
	np.random.seed()
	
	localUserTopicPref = list(np.random.dirichlet(hyperGamma))
	userTopicPrefVectors[x] = localUserTopicPref

	writeStr = str(x) + ' ' + ' '.join([str(i) for i in localUserTopicPref]) 

	userTopicPrefFile.write(writeStr + "\n")

userTopicPrefFile.close()
'''
print("Got user topic preferences", len(userTopicPrefVectors))


# print("Read data from file")
print("Generating events...")
allSyntheticEvents = generate_synthetic_events()
print(len(allSyntheticEvents))

writeOnlyEventsToFile()

allSyntheticDocs = generate_synthetic_docs(allSyntDocsFileName)

print "Generated ", len(allSyntheticEvents), " events and docs....\n"