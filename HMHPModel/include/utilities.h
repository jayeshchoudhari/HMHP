#ifndef UTILITIES_H 
#define UTILITIES_H

int getSampleFromMultinomial(std::vector<double> calculatedProbVec);
std::vector<double> getNormalizedLogProb(std::vector<double> calculatedProbVec);
int getSampleFromDiscreteDist(std::vector<double> normalizedProbVector);
double getSampleFromGamma(double alpha, double beta);

#endif