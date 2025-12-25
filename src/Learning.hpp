#ifndef LEARNING_HPP
#define LEARNING_HPP

#include "Types.hpp"
#include <unordered_map>
#include <vector>

std::vector<double> learnability(
    const Repertoire& repertoire,
    const Parents& parents,
    const AdjacencyMatrix& adjMatrix
);


double stayProbability(std::vector<std::pair<Repertoire, double>> transitions);

std::vector<double> baseWeights(
    Strategy strategy,
    const Repertoire& repertoire,
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope,
    const std::vector<double>& statePayoffs   
);

std::pair<std::vector<Repertoire>, std::vector<std::vector<std::pair<Repertoire, double>>>>  generateReachableRepertoires(
    Strategy strategy, 
    const AdjacencyMatrix& adjMatrix, 
    const PayoffVector& payoffs, 
    const std::vector<double>& traitFrequencies,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope,
    const std::vector<double>& statePayoffs,
    double transparency,
    const std::vector<double>& depths,
    TraitDistribution dist
);

std::vector<Repertoire> generateAllRepertoires(const AdjacencyMatrix& adjMatrix);

size_t countLearnedTraits(const Repertoire& r);

std::vector<Repertoire> retrieveBetterRepertoires(const std::vector<Repertoire>& repertoires, const Repertoire& singleRepertoire);

#endif // LEARNING_HPP