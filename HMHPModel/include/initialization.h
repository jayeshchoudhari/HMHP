#ifndef INITIALIZATION_H 
#define INITIALIZATION_H

class InitializeModel
{
	private:

	public:

		double baseAlpha;
		double baseBeta; 
		double logBase;

		InitializeModel(DataIO &dataIOObj);
		int initializeForSampler(DataIO &dataIOObj);
		int updateCountMatrices(DataIO &dataIOObj, std::vector< std::vector <li> > &localNewSyntheticEvents);
		int updateNodeNodeCountMap(DataIO &dataIOObj);
		int initializeUserUserInfluence(DataIO &dataIOObj);
		int initializeAvgProbabilityVectors(DataIO &dataIOObj);
		int initializeAvgTopicProbabilityVectors(DataIO &dataIOObj);
};


#endif