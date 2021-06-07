from __future__ import division
from collections import defaultdict
from collections import Counter
import sys
import math
import numpy as np
# import matplotlib.pyplot as plt


class SyntheticDataGenerator(object):
    def __init__(self, tMax, seedVal, eventsToBeGenerated,
                 numTopics, alphaMultiplier, maxLevelsToBeGenerated = 4):
        """
        Init Method
        """
        # print("inside init -- ", tMax)
        self.tMax = tMax,
        self.seedVal = seedVal
        self.eventsToBeGenerated = eventsToBeGenerated
        self.numTopics = numTopics
        self.alphaMultiplier = alphaMultiplier
        self.maxLevelsToBeGenerated = maxLevelsToBeGenerated

    def get_homogenous_pp_timestamps(self, lambda_rate):
        """
        DESC
        """
        hppTimeStamps = []
        prev_ti = 0
        # print("inside get-homoge-pp-tmax - ", self.tMax)
        # print("inside get-homoge-pp- seed -- ", self.seedVal)
        max_Time = self.tMax[0]
        while prev_ti <= max_Time:
                key_u = np.random.uniform(0,1)
                key_x = prev_ti - (np.log(key_u)/lambda_rate)
                if key_x <= max_Time:
                        hppTimeStamps.append(key_x)
                prev_ti = key_x
        return hppTimeStamps

    def get_lambda_for_level_l(self, startTime, x, W_uv):
        """
        DESC
        """
        timeDiff = x - startTime
        # func_value = W_uv * math.exp(-timeDiff)
        func_value = W_uv * math.exp(- self.alphaMultiplier * timeDiff)
        return func_value

    def get_nhpp_timestamps_for_level_l(self, lambdaUpperBound, t_e, W_uv):
        """
        DESC
        """
        T, t0, n, m = self.tMax, t_e, 0, 0
        lambda_constant = lambdaUpperBound
        sm = t_e
        # homoTime = []
        inhomoTime = []
        # allDs = []
        while sm < T:
            u = np.random.uniform(0,1)
            w = -(np.log(u) / lambda_constant)  # so that w ~ exp(lambda_constant)
            sm = sm + w
            probVal = (self.get_lambda_for_level_l(t_e, sm, W_uv))/lambda_constant
            # probVal = (get_lambda_for_level_l_rayleigh(t_e, sm, W_uv))/lambda_constant
            if sm < T and probVal >= 1e-7:
                # homoTime.append(sm)
                # sm are the points in the homo. PP
                d = np.random.uniform(0,1)
                # print(sm, d, math.exp(-(sm - t_e)))
                if d <= probVal:
                    # allDs.append(d)
                    inhomoTime.append(sm)   # inhomoTime are the points in inhomo. PP
                    # n += 1
                    # m += 1
                else:
                    break
        # print (inhomoTime)
        # plot_for_confirmation(homoTime, inhomoTime, allDs)
        return inhomoTime

    def generate_synthetic_events(self, userBaseRates, userTopicPrefVectors,
            userInfluence, topicTopicProbVectors):
        """
        Processing Method
        @userBaseRates (list): desc
        @userTopicPrefVectors (list): desc
        @userInfluence (list): desc
        @topicTopicProbVectors (list): desc
        return: Desc
        """
        print("Tmax -- ", self.tMax)
        print("Seed -- ", self.seedVal)

        np.random.seed(self.seedVal)
        topics = [i for i in range(0, self.numTopics)]
        W = userInfluence
        # numNodes = len(nodes)
        noparent = '-1'
        allEvents = []
        # will be a sorted list of events
        # where each event is a list [time, node, parent, topic]
        # topic to each doc/event will be appended later...
        eventsWithId = {}
        level0Events = []
        # Generate Level-0 events
        for node in userBaseRates.keys():
            # utils.get_nhpp_timestamps(Tmax, lambdaFunction, lambdaUpperBound) 
            # where lambdaUpperBound is a constant lambda upper bounding the rate
            # or intensity function
            # currently we are assuming that is known for each node                  

            userMuVal = userBaseRates[node]
            nodeL0Timestamps = self.get_homogenous_pp_timestamps(userMuVal)
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

        eventsCount = len(allEvents)
        # Generate level-l events
        prevLevelEvents = level0Events
        currentLevelEvents = []
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
                # if node u is present in the followers map or if u has followers
                # then go for generating events from the neighbors of u...
                # if type(upresent) == dict:
                if upresent != -1:
                    for node in W[u]:
                        W_uv = W[u][node]
                        # nhppTimestampsLevelL = utils.get_nhpp_timestamps_for_level_l(Tmax,
                        #                        1, u_ts, W_uv)
                        # nhppTimestampsLevelL = utils.get_nhpp_timestamps_for_level_l(Tmax,
                        #                         lambda_constant, u_ts, W_uv)
                        if W_uv > 0:
                            nhppTimestampsLevelL = self.get_nhpp_timestamps_for_level_l(W_uv,
                                                                            u_ts, W_uv)
                            srcDest = str(u) + "_" + str(node)
                            inteArrivSum = 0
                            for ts in nhppTimestampsLevelL:
                                probVec = np.array(topicTopicProbVectors[parentTopic])
                                if probVec.sum() == 0:
                                    eta_n = np.random.choice(topics, 1)[0]
                                else:
                                    probVec /= probVec.sum().astype(float)
                                    eta_n = np.random.choice(topics, 1, True, probVec)[0]
                                currentLevelEvents.append([ts, node, parentEventId,
                                                           eta_n, level])
                                eventsCount += 1
                        if eventsCount > self.eventsToBeGenerated:
                            nomoreEventsFlag = 0
                            break
                if nomoreEventsFlag == 0 :
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
            if level > self.maxLevelsToBeGenerated:
                print("done with required levels... Breaking here....")
                break
            if len(allEvents) > self.eventsToBeGenerated:
                print("Events more than required.. Breaking here")
                break
        allEvents.sort(key = lambda row:row[0])
        uniqEventIdSortedIndex = dict()
        for eid in range(len(allEvents)):
            ts = allEvents[eid][0]
            node = allEvents[eid][1]
            uniqEventId = str(node) + "_" + str(ts)
            uniqEventIdSortedIndex[uniqEventId] = eid
        # print("Adding the required parent Index to the events...")
        # now for each event, if its parent is some uniqEventId,
        # we have change that with the event index....
        for eid in range(len(allEvents)):
            uniqParentId = allEvents[eid][2]
            if uniqParentId != '-1':
                indexOfParent = uniqEventIdSortedIndex[uniqParentId]
                allEvents[eid][2] = indexOfParent
            else:
                allEvents[eid][2] = -1
        print("numEvents = ", len(allEvents))
        return allEvents