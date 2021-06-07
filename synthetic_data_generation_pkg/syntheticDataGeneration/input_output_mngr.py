import os
import sys
import numpy as np

class FileInputManager(object):

    def __check_file_exist(self, file_path):
        if not os.path.exists(file_path):
            raise Exception("File Not Found. PATH: {}"
                            .format(file_path))

    def getUserBaseRates(self, userBaseRatesFilePath):
        """
        Method to read the file with path userBaseRatesFilePath
        @userBaseRatesFilePath (str): User Base rate file path
        return: (dict)
        """
        tempUBR = dict()
        self.__check_file_exist(userBaseRatesFilePath)
        with open(userBaseRatesFilePath, 'r') as fp:
            for line in fp:
                splts = line.strip().split()
                uid = int(splts[0].strip())
                ubr = float(splts[1].strip())
                tempUBR[uid] = ubr
        return tempUBR

    def getUserUserInfluence(self, userUserInfluenceFilePath):
        """
        Method to read the file with path userUserInfluenceFilePath and return dict
        @userUserInfluenceFilePath (str): User Influence File Path
        return: (dict)
        """
        tempWuv = dict()
        self.__check_file_exist(userUserInfluenceFilePath)
        with open(userUserInfluenceFilePath, 'r') as wuvfile:
            for line in wuvfile:
                splitLine = line.split()
                # count = int(splitLine[0].strip())
                uNode = int(splitLine[1].strip())
                tempWuv[uNode] = {}
                j = 2
                while j < len(splitLine):
                    vNode = int(splitLine[j].strip())
                    uvInf = float(splitLine[j+1].strip())
                    tempWuv[uNode][vNode] = uvInf
                    j = j + 2
        return tempWuv

    def getUserTopicPrefVectors(self, userTopicPrefVectorsFilePath):
        """
        Method to read the file with path userTopicPrefVectorsFilePath and return dict
        @userTopicPrefVectorsFilePath (str): User Topic Preferences File Path
        return: (dict)
        """
        tempUserTopicPrefVector = dict()
        self.__check_file_exist(userTopicPrefVectorsFilePath)
        with open(userTopicPrefVectorsFilePath, 'r') as fp:
            for line in fp:
                splts = line.strip().split()
                uid = int(splts[0].strip())
                tempVec = []
                # sumVec = 0
                for j in range(1,len(splts)):
                    ele = float(splts[j].strip())
                    # sumVec += ele
                    tempVec.append(ele)
                tempVec = np.array(tempVec)
                tempVec /= tempVec.sum()
                tempUserTopicPrefVector[uid] = list(tempVec)
        return tempUserTopicPrefVector

    def getTopicTopicProbVectors(self, topicTopicProbVectorsFilePath):
        """
        Method to read the file with path topicTopicProbVectorsFilePath and return list of lists
        @topicTopicProbVectorsFilePath (str): Topic Topic influence/interaction File Path
        return: list of lists
        """
        tempTopicTopicVec = []
        self.__check_file_exist(topicTopicProbVectorsFilePath)
        with open(topicTopicProbVectorsFilePath, 'r') as fp:
            for line in fp:
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
        return tempTopicTopicVec

    def getTopicWordProbVectors(self, wordDistTopicsFilePath):
        """
        Method to read the file with path wordDistTopicsFilePath and return list of lists
        @wordDistTopicsFilePath (str): Topic word distribution File Path
        return: list of lists
        """
        tempWordDistVectors = []
        self.__check_file_exist(wordDistTopicsFilePath)
        with open(wordDistTopicsFilePath, 'r') as fp:
            for line in fp:
                splts = line.strip().split()
                # tid = int(splts[0].strip())
                tempVec = []
                for j in range(1, len(splts)):
                    tempVec.append(float(splts[j].strip()))
                tempVec = np.array(tempVec)
                tempVec /= tempVec.sum()
                tempWordDistVectors.append(list(tempVec))
        return tempWordDistVectors



class FileOutputManager(object):
    def __init__(self, allSyntheticEvents):
        self.allSyntheticEvents = allSyntheticEvents

    def writeOnlyEventsToFile(self, allSyntEventsFilePath):
        """
        Method to write all the events to the file at allSyntEventsFilePath
        @allSyntEventsFilePath (str): write events to the output file 
        return: None
        """
        with open(allSyntEventsFilePath, 'w') as fp:
            # Sample topic for each event based on the parent
            for event in range(0, len(self.allSyntheticEvents)):
                eventStr = ' '.join([str(ele) for ele in self.allSyntheticEvents[event]])
                eventStr = eventStr + "\n"
                fp.write(eventStr)
                fp.flush()

    def generate_synthetic_docs(self, vocabsize, wordDistTopics, allSyntDocsFilePath):
        """
        Method to generate documents write all the documents corresponding to the events to the file at allSyntDocsFilePath
        @vocabsize (int): Describes the maximum number of words in the vocabulary
        @wordDistTopics (list of lists): Describes topic word distribtion or distribution over vocabulary for each topic 
        @allSyntDocsFilePath (str): write documents to the output file 
        return: None
        """
        print("Generating Documents...")
        vocabulary = [i for i in range(0, vocabsize)]
        # Generate word distribution for each topic
        # hyperAlpha = [np.random.uniform(0,1) for k in range(0,vocabsize)]
        # Generate words for all the docs and write to file
        docWords = []
        print(allSyntDocsFilePath)
        with open(allSyntDocsFilePath, 'w') as fp:
            # for eachdoc in allSyntheticEventsTopicsAssigned:
            for eachdoc in self.allSyntheticEvents:
                docsize = np.random.poisson(7)
                wordsInDoc = np.random.choice(vocabulary,
                                              docsize,
                                              True,
                                              wordDistTopics[eachdoc[3]])
                # draw doc size number of words
                docWords.append(wordsInDoc)
                docStr = ' '.join([str(word) for word in wordsInDoc])
                docStr = docStr + "\n"
                fp.write(docStr)
                fp.flush()
        return docWords
