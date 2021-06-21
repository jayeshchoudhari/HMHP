#ifndef DATAIO_H 
#define DATAIO_H

class DataIO
{
    private:

    public:
        // global members declaration
        vector< vector <li> > allEvents;
        vector< vector <li> > newSyntheticEvents;
        vector< vector <ui> > allDocs;
        vector < vector <ui> > wordHistAllDocsVector;
        vector < vector <ui> > followersMap;
        vector < vector <ui> > reverseFollowersMap;

        vector <vector <ui> > allPossibleParentEvents;
        vector <vector <double> > allPossibleParentEventsExponentials;

        string inputFilePathsFileName, outputFilePathsFileName;
        unordered_map <string, string> configInputFiles, configOutputFiles;

        DataIO(string ipFilePathsFileName, string opFilePathsFileName);
        unordered_map <string, string> getConfigInputOutputFileNames(string fileName);
        int populateDataInDataStructure(unordered_map<string, string> configInputFiles);
        
        vector <vector <ui> > populateParentEventsIds(string parentEventsFilename);
        vector <vector <double> > populateParentEventsExpKernelValue(string parentEventExpKernelValFilename);
        
        vector < vector <ui> > getSyntheticDocsFromFile(string fileName);
        int DataIO :: createWordHistForAllDocs();
};


#endif