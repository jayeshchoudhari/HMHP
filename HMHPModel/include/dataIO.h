#ifndef DATAIO_H 
#define DATAIO_H

class DataIO
{
    private:

    public:
        // global members declaration
        std::string inputFilePathsFileName, outputFilePathsFileName;
        std::unordered_map <std::string, std::string> configInputFiles, configOutputFiles;

        ui maxLevel, totalWords, maxNumNodes, numTopics, vocabSize;

        std::vector< std::vector <li> > allEvents;
        std::vector< std::vector <li> > newSyntheticEvents;
        
        std::vector< std::vector <ui> > allDocs;
        std::vector < std::vector <ui> > wordHistAllDocsVector;
        
        std::vector < std::vector <ui> > followersMap;
        std::vector < std::vector <ui> > reverseFollowersMap;

        std::vector <std::vector <ui> > allPossibleParentEvents;
        std::vector <std::vector <double> > allPossibleParentEventsExponentials;

        std::map <ui, std::vector<ui> > childEventsMap;
        std::map <ui, std::vector<ui> > nodeEventsMap;
        std::vector<ui> nodeEventsCountMap;

        std::vector<ui> levelInfo;


        std::unordered_map<int, int> outDegreeMap, inDegreeMap;

        std::unordered_map<ui, double> eventIndexTimestamps;
        std::unordered_map <ui, int> eventIndexDocSizeMap;


        std::vector<double> hyperBeta, hyperGamma, hyperAlpha;
        double sumAlpha, sumBeta, sumGamma;

        std::map<ui, std::map<ui, ui> > nodeNodeCount;
        std::unordered_map<ui, std::vector <ui> > nodeNodeCombUpdateInfluence;

        std::vector < std::vector < ui > > NTWCountVec;
        std::vector < std::vector < ui > > NTTCountVec;
        std::vector < std::vector < ui > > NUTCountVec;

        std::vector < ui > NTWSumWordsVec;
        std::vector < ui > NTTSumTopicsVec;
        std::vector < ui > NUTSumTopicsVec;

        std::unordered_map <ui, std::unordered_map<ui, double> > userUserInfluence;

        std::unordered_map<std::string, int> groupTransactionsSum, groupSourceTransactionsSum, actualEdgesNuSum;

        std::vector <std::vector <double> > avgProbParForAllEvents;
        std::vector <std::vector <double> > avgTopicProbVector;


        DataIO(std::string ipFilePathsFileName, std::string opFilePathsFileName, int numOfTopics, int maxNumOfNodes);
        std::unordered_map <std::string, std::string> getConfigInputOutputFileNames(std::string fileName);
        int populateDataInDataStructure(std::unordered_map <std::string, std::string> configInputFiles);
        
        std::vector < std::vector <li> > getSyntheticEventsFromFile(std::string fileName);
        std::vector < std::vector <ui> > getSyntheticDocsFromFile(std::string fileName);

        std::vector <std::vector <ui> > populateParentEventsIds(std::string parentEventsFilename);
        std::vector <std::vector <double> > populateParentEventsExpKernelValue(std::string parentEventExpKernelValFilename);
        
        std::unordered_map<int, int> getDegreeMap(std::string degreeMapFile);
        int readHyperParametersVectorsFromFile();

        int createWordHistForAllDocs();


        std::vector < std::vector <ui> > readIntVectorMapFromFile(std::string fileName);

        std::unordered_map <ui, ui> getHistOfWordsOverWordsFromDoc(std::vector<ui> doc);
};


#endif