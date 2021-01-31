"""
Synthetic Data Generation Tool

Author: Jayesh Choudhari <choudhari.jayesh@iitgc.ac.in>
Date: 30 Jan, 2021
"""



import argparse
import os

from syntheticDataGeneration.input_output_mngr import (FileInputManager,
                                                       FileOutputManager)
from syntheticDataGeneration.data_generator import SyntheticDataGenerator

def argument_parser():
    """
    Method to create Project Argumeht parser
    """
    parser = argparse.ArgumentParser(
            description='Synthetic Data Generation Tool')
    parser.add_argument('--tmax', '-t', type=int,
                         required=True,
                         metavar='\b',
                         dest='tMax',
                         help='Tmax (time horizon)')
    parser.add_argument('--random-seed', '-s', type=int,
                        required=True,
                        metavar='\b',
                        dest='seedVal',
                        help='random_seed value')
    parser.add_argument('--output-event-count', '-e', type=int,
                        metavar='\b',
                        required=True,
                        dest='eventsToBeGenerated',
                        help='max number of events to be generated')
    parser.add_argument('--topic-count', '-c', type=int,
                        metavar='\b',
                        dest='numTopics',
                        help='no. of topics')
    parser.add_argument('--vocab-size', '-v', type=int,
                        metavar='\b',
                        dest='vocabSize',
                        help='Vocab Size')
    parser.add_argument('--alpha-multiplier', '-a', type=int,
                        metavar='\b',
                        dest='alphaMultiplier',
                        default=1,
                        help='alphaMultiplier. This variable is optional. '
                             'currently not used')
    parser.add_argument('--inputDir', '-i', type=str,
                        metavar='\b',
                        dest='inputDir',
                        default= os.path.dirname(os.getcwd()) + '/inputData/',
                        help='Input Directory path.\n' +
                             'Drirectory  must present below file. '
                             '1. sampleUsersBaseRate.txt '
                             '2. sampleUserUserInfluence.txt '
                             '3. topicTopicDistribution.txt '
                             '4. topicWordDistribution.txt '
                             '5. userTopicDistribution.txt')
    parser.add_argument('--otputDir', '-o', type=str,
                        metavar='\b',
                        dest='outputDir',
                        default= os.path.dirname(os.getcwd()) + '/outputData/',
                        help='Output Directory path.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = argument_parser()
    input_obj = FileInputManager()

    INPUT_DIR_PATH = args.inputDir
    if not INPUT_DIR_PATH.endswith('/'):
        INPUT_DIR_PATH += '/'

    OUTPUT_DIR_PATH = args.outputDir
    if not OUTPUT_DIR_PATH.endswith('/'):
        OUTPUT_DIR_PATH += '/'

    userBaseRatesFileName = 'sampleUsersBaseRate.txt'
    userBaseRatesFilePath = INPUT_DIR_PATH + userBaseRatesFileName

    userUserInfluenceFileName = 'sampleUserUserInfluence.txt'
    userUserInfluenceFilePath = INPUT_DIR_PATH + userUserInfluenceFileName

    userTopicPrefVectorsFileName = 'userTopicDistribution.txt'
    uuserTopicPrefVectorsFilePath = INPUT_DIR_PATH + userTopicPrefVectorsFileName

    topicTopicProbVectorsFileName = 'topicTopicDistribution.txt'
    topicTopicProbVectorsFilePath = INPUT_DIR_PATH + topicTopicProbVectorsFileName

    wordDistTopicsFileName = 'topicWordDistribution.txt'
    wordDistTopicsFilePath = INPUT_DIR_PATH + wordDistTopicsFileName
    #1. ------------------------------------------------------------------------------
    print("Reading Input Files...")
    user_base_rates =  input_obj.getUserBaseRates(userUserInfluenceFilePath)
    user_topic_vtr = input_obj.getUserTopicPrefVectors(userTopicPrefVectorsFilePath)
    user_influence = input_obj.getUserUserInfluence(userUserInfluenceFilePath)
    topic_pro_vtr = input_obj.getTopicTopicProbVectors(topicTopicProbVectorsFilePath)

    #2. ------------------------------------------------------------------------------
    print("Generating Data...")
    sdc_obj = SyntheticDataGenerator(tMax=arg.tMax,
                                     seedVal=args.seedVal,
                                     eventsToBeGenerated=args.eventsToBeGenerated,
                                     numTopics=args.numTopics,
                                     alphaMultiplier=args.alphaMultiplier)
    allSyntheticEvents = sdc_obj.generate_synthetic_events(userBaseRates=user_base_rates,
                                                    userTopicPrefVectors=user_topic_vtr,
                                                    userInfluence=user_influence,
                                                    topicTopicProbVectors=topic_pro_vtr)

    #3. ------------------------------------------------------------------------------
    print("Writing Data...")
    allSyntEventsFileName = 'eventsFile.txt'
    allSyntEventsFilePath = OUTPUT_DIR_PATH + allSyntDocsFileName

    allSyntDocsFileName = 'documentsFile.txt'
    allSyntDocsFilePath = OUTPUT_DIR_PATH + allSyntDocsFileName
    wordDistTopics = input_obj.getTopicWordProbVectors(wordDistTopicsFilePath)

    output_obj = FileOutputManager(allSyntheticEvents)
    output_obj.writeOnlyEventsToFile(allSyntEventsFilePath)
    resp = output_obj.generate_synthetic_docs(vocabsize=args.vocabSize,
                                              wordDistTopics=wordDistTopics,
                                              allSyntDocsFilePath=allSyntEventsFilePath)

