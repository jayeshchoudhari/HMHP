#include "namespace.h"
#include "utilities.h"

using namespace std;

/*
// This is prone to overflow....
vector<double> getNormalizedLogProb(vector<double> calculatedProbVec)
{
	vector<double> normalizedProbVector(calculatedProbVec.size(), 0.0);

	// each of the terms will be normalized as:
	// log(pi) - ( log(mi) + log(sum (exp (log(pi) - log(mi)))) )

	// double maxTerm = *max_element(calculatedProbVec.begin(), calculatedProbVec.end());

	double sumOverExp = 0.0;
	// sum(exp(log(pi) - log(mi)))
	for(unsigned int i = 0; i < calculatedProbVec.size(); i++)
	{
		// sumOverExp += exp(calculatedProbVec[i] - maxTerm);
		sumOverExp += exp(calculatedProbVec[i]);
	}

	for(unsigned int i = 0; i < calculatedProbVec.size(); i++)
	{
		// normalizedProbVector[i] = calculatedProbVec[i] - (log(maxTerm) + log(sumOverExp));
		// normalizedProbVector[i] = calculatedProbVec[i] - (log(sumOverExp));
		normalizedProbVector[i] = exp(calculatedProbVec[i] - (log(sumOverExp)));
	}

	return normalizedProbVector;
}
*/

vector<double> getNormalizedLogProb(vector<double> calculatedProbVec)
{
	vector<double> normalizedProbVector(calculatedProbVec.size(), 0.0);

	// each of the terms will be normalized as:
	// log(pi) - ( log(mi) + log(sum (exp (log(pi) - log(mi)))) )

	double maxTerm = *max_element(calculatedProbVec.begin(), calculatedProbVec.end());

	double sumOverExp = 0.0;
	// sum(exp(log(pi) - log(mi)))
	for(unsigned int i = 0; i < calculatedProbVec.size(); i++)
	{
		sumOverExp += exp(calculatedProbVec[i] - maxTerm);
	}

	for(unsigned int i = 0; i < calculatedProbVec.size(); i++)
	{
		normalizedProbVector[i] = exp(calculatedProbVec[i] - (maxTerm + log(sumOverExp)));
		// normalizedProbVector[i] = exp(calculatedProbVec[i] - (log(sumOverExp)));
	}

	return normalizedProbVector;
}

int getSampleFromDiscreteDist(vector<double> normalizedProbVector)
{
	// remove .. later
	random_device rd;
	mt19937 gen(rd());
	
	// default_random_engine gen;

	discrete_distribution<int> distribution(normalizedProbVector.begin(), normalizedProbVector.end());
	int ind = distribution(gen);

	return ind;
}


double getSampleFromGamma(double alpha, double beta)
{
	random_device rd;
	mt19937 gen(rd());
	gamma_distribution<double> distribution(alpha, beta);

	double gammaVal = distribution(gen);

	return gammaVal;
}
