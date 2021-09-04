#ifndef GIBBSSAMPLER_H 
#define GIBBSSAMPLER_H

class GibbsSampler
{
	private:

	public:

		double baseAlpha;
		double baseBeta; 
		double logBase;
		double logLikelihood;


		// Sampler Initialization
		GibbsSampler(DataIO &dataIOObj);
		int initializeForSampler(DataIO &dataIOObj);
		int updateCountMatrices(DataIO &dataIOObj, std::vector< std::vector <li> > &localNewSyntheticEvents);
		int updateNodeNodeCountMap(DataIO &dataIOObj);
		int initializeUserUserInfluence(DataIO &dataIOObj);
		int initializeAvgProbabilityVectors(DataIO &dataIOObj);
		int initializeAvgTopicProbabilityVectors(DataIO &dataIOObj);


		// Iterative MCMC Sampler
		int IterativeSampler(DataIO &dataIOObj, ui ITERATIONS, int BURN_IN);


		// Topic Inference
		int sampleTopicAssignment(DataIO &dataIOObj, int ITE, int BURN_IN);
		li getSampledTopicAssignment(DataIO &dataIOObj, li eventIndex, li eventNode, li eventParent, std::vector<ui> doc, int ITE, int BURN_IN);
		double getFirstTermOfTopicAssignmentCondProb(DataIO &dataIOObj, li eventNode, li eventParent, li topic);
		double getMiddleTermOfTopicAssignmentCondProb(DataIO &dataIOObj, li topic, std::unordered_map <int, ui> childEventTopicsHist, li eventIndex, li eventParent);
		std::unordered_map <int, ui>  getHistOfTopicsOverChildEvents(DataIO &dataIOObj, li eventIndex);
		double  getThirdTermOfTopicAssignmentCondProb(DataIO &dataIOObj, li topic, std::vector<ui> doc, std::vector<ui> wordHistVec);
		int  printTopicTopicCount(DataIO &dataIOObj);


		// Parent Inference
		int sampleParentAssignment(DataIO &dataIOObj, int ITE, int BURN_IN);
		li getSampledParentAssignment(DataIO &dataIOObj, double eventTime, li eventNode, li eventIndex, li eventTopic, int ITE, int BURN_IN);
		std::vector<double> populateCalculatedProbVec(DataIO &dataIOObj, std::vector <ui> possibleParentEvents, std::vector<double> possibleParentExp, li eventNode, li eventTopic, double eventTime, int ITE);
		double getFirstTermOfParentAssignment(DataIO &dataIOObj, li possParentEventTopic, li eventTopic);
		double getFirstTermOfParentAssignmentNoParent(DataIO &dataIOObj, li eventNode, li eventTopic);
		std::vector<ui> getPossibleParentEvents(DataIO &dataIOObj, li eventNode, li eventIndex);


		// User-User Influence Inference
		int sampleInfluenceAssignment(DataIO &dataIOObj);


		// UserBaseRate 
		int updateUserBaseRates(DataIO &dataIOObj);


		// Count Matrices Handler...

		// Decrement Counters...
		int decreamentCountFromMatrices(DataIO &dataIOObj, li eventIndex, li eventNode, li eventParent, li eventTopic, std::vector <ui> doc, bool topicSampling);
		int decreamentCountFromTopicTopic(DataIO &dataIOObj, li eventParentTopic, li eventTopic);
		int decreamentCountFromUserTopic(DataIO &dataIOObj, li eventNode, li eventTopic);
		int decreamentCountFromTopicWord(DataIO &dataIOObj, li eventTopic, std::vector <ui> doc);
		int decreamentCountsFromChildEvents(DataIO &dataIOObj, li eventTopic, li eventIndex);

		// Increment Counters...
		int increamentCountToMatrices(DataIO &dataIOObj, li eventIndex, li eventNode, li eventParent, li eventTopic, std::vector<ui> doc, bool topicSampling);
		int increamentCountInTopicTopic(DataIO &dataIOObj, li eventParentTopic, li eventTopic);
		int increamentCountInUserTopic(DataIO &dataIOObj, li eventNode, li eventTopic);
		int increamentCountInTopicWord(DataIO &dataIOObj, li eventTopic, std::vector<ui> doc);
		int increamentCountsForChildEvents(DataIO &dataIOObj, li eventTopic, li eventIndex);


		int getSampleFromMultinomial(std::vector<double> calculatedProbVec);

};


#endif