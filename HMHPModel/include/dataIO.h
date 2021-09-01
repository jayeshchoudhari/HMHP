#ifndef DATAIO_H 
#define DATAIO_H

class DataIO
{
    private:

    public:
        // global members declaration
        std::vector< std::vector <li> > allEvents;
        std::vector< std::vector <li> > newSyntheticEvents;
        std::vector< std::vector <ui> > allDocs;
        std::vector < std::vector <ui> > wordHistAllDocsVector;
        std::vector < std::vector <ui> > followersMap;
        std::vector < std::vector <ui> > reverseFollowersMap;

        std::vector <std::vector <ui> > allPossibleParentEvents;
        std::vector <std::vector <double> > allPossibleParentEventsExponentials;

        std::string inputFilePathsFileName, outputFilePathsFileName;
        std::unordered_map <std::string, std::string> configInputFiles, configOutputFiles;

        std::unordered_map<int, int> outDegreeMap, inDegreeMap;

        std::unordered_map<ui, double> eventIndexTimestamps;
        std::unordered_map <ui, int> eventIndexDocSizeMap;

        std::vector<ui> levelInfo;

        ui maxLevel, totalWords, maxNumNodes, numTopics, vocabSize;

        std::vector<double> hyperBeta, hyperGamma, hyperAlpha;
        double sumAlpha, sumBeta, sumGamma;


        DataIO(std::string ipFilePathsFileName, std::string opFilePathsFileName);
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