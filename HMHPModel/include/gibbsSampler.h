#ifndef GIBBSSAMPLER_H 
#define GIBBSSAMPLER_H

class GibbsSampler
{
	private:

	public:

		double baseAlpha;
		double baseBeta; 
		double logBase;

		GibbsSampler(DataIO &dataIOObj);
		int initializeForSampler(DataIO &dataIOObj);
		int updateCountMatrices(DataIO &dataIOObj, std::vector< std::vector <li> > &localNewSyntheticEvents);
		int updateNodeNodeCountMap(DataIO &dataIOObj);
		int initializeUserUserInfluence(DataIO &dataIOObj);
		int initializeAvgProbabilityVectors(DataIO &dataIOObj);
		int initializeAvgTopicProbabilityVectors(DataIO &dataIOObj);

		int sampleInfluenceAssignment(DataIO &dataIOObj);
		int updateUserBaseRates(DataIO &dataIOObj);

};


#endif