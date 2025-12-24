#ifndef LEARNING_HPP
#define LEARNING_HPP

#include "Types.hpp"
#include <unordered_map>
#include <vector>
#include <string>

std::vector<double> learnability(
    const Repertoire& repertoire,
    const Parents& parents,
    const AdjacencyMatrix& adjMatrix
);


double stayProbability(std::vector<std::pair<Repertoire, double>> transitions);

struct RepertoireHash {
    std::size_t operator()(const Repertoire& repertoire) const {
        return std::hash<std::string>{}(std::string(repertoire.begin(), repertoire.end()));
    }
};

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
    double transparency
);

std::vector<Repertoire> generateAllRepertoires(const AdjacencyMatrix& adjMatrix);

size_t countLearnedTraits(const Repertoire& r);

std::vector<Repertoire> retrieveBetterRepertoires(const std::vector<Repertoire>& repertoires, const Repertoire& singleRepertoire);

#endif // LEARNING_HPP