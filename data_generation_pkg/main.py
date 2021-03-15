"""
Synthetic Data Generation Tool
Author: Jayesh Choudhari <choudhari.jayesh@alumni.iitgc.ac.in>
"""



import argparse
import os

from syntheticDataGeneration.input_output_mngr import (FileInputManager,
                                                       FileOutputManager)
from syntheticDataGeneration.data_generator import SyntheticDataGenerator

def argument_parser():
    """
    Method to create Project Argument parser
    """
    parser = argparse.ArgumentParser(
            description='Synthetic Data Generation Tool')
    parser.add_argument('--tmax', '-t', type=int,
                         required=True,
                         metavar='\b',
                         dest='tMax',
                         help='Tmax (time horizon)')
    parser.add_argument('--seed', '-s', type=int,
                        required=True,
                        metavar='\b',
                        dest='seedVal',
                        help='random_seed value -- fix this value to have the same output every time...')
    parser.add_argument('--output-event-count', '-e', type=int,
                        required=True,
                        metavar='\b',
                        dest='eventsToBeGenerated',
                        help='max number of events to be generated')
    parser.add_argument('--topic-count', '-c', type=int,
                        required=True,
                        metavar='\b',
                        dest='numTopics',
                        help='max no. of topics')
    parser.add_argument('--vocab-size', '-v', type=int,
                        required=True,
                        metavar='\b',
                        dest='vocabSize',
                        help='Vocab Size')
    parser.add_argument('--alpha-multiplier', '-a', type=int,
                        metavar='\b',
                        dest='alphaMultiplier',
                        default=1,
                        help='alphaMultiplier. For now, this variable is optional. '
                             'currently not used')
    parser.add_argument('--inputDir', '-i', type=str,
                        metavar='\b',
                        dest='inputDir',
                        default= os.getcwd() + '/inputData/',
                        help='Input Directory path.\n' +
                             'Directory must consist below files. '
                             '1. sampleUsersBaseRate.txt '
                             '2. sampleUserUserInfluence.txt '
                             '3. topicTopicDistribution.txt '
                             '4. topicWordDistribution.txt '
                             '5. userTopicDistribution.txt')
    parser.add_argument('--otputDir', '-o', type=str,
                        metavar='\b',
                        dest='outputDir',
                        default= os.getcwd() +  '/outputData/',
                        help='Output Directory path.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    print(os.path.dirname(os.getcwd()))
    args = argument_parser()
    input_obj = FileInputManager()

    INPUT_DIR_PATH = args.inputDir
    print(INPUT_DIR_PATH)
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
    userTopicPrefVectorsFilePath = INPUT_DIR_PATH + userTopicPrefVectorsFileName

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
    print("Printing all args --- ")
    print(args.tMax)
    sdg_obj = SyntheticDataGenerator(tMax=args.tMax,
                                     seedVal=args.seedVal,
                                     eventsToBeGenerated=args.eventsToBeGenerated,
                                     numTopics=args.numTopics,
                                     alphaMultiplier=args.alphaMultiplier)
    allSyntheticEvents = sdg_obj.generate_synthetic_events(userBaseRates=user_base_rates,
                                                    userTopicPrefVectors=user_topic_vtr,
                                                    userInfluence=user_influence,
                                                    topicTopicProbVectors=topic_pro_vtr)

    #3. ------------------------------------------------------------------------------
    print("Writing Data...")
    allSyntEventsFileName = 'eventsFile.txt'
    allSyntEventsFilePath = OUTPUT_DIR_PATH + allSyntEventsFileName

    allSyntDocsFileName = 'documentsFile.txt'
    allSyntDocsFilePath = OUTPUT_DIR_PATH + allSyntDocsFileName
    wordDistTopics = input_obj.getTopicWordProbVectors(wordDistTopicsFilePath)

    output_obj = FileOutputManager(allSyntheticEvents)
    output_obj.writeOnlyEventsToFile(allSyntEventsFilePath)
    resp = output_obj.generate_synthetic_docs(vocabsize=args.vocabSize,
                                              wordDistTopics=wordDistTopics,
                                              allSyntDocsFilePath=allSyntDocsFilePath)
