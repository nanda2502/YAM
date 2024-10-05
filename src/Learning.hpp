#ifndef LEARNING_HPP
#define LEARNING_HPP

#include "Types.hpp"
#include <vector>
#include <string>
#include <random>

std::vector<bool> learnability(const Repertoire& repertoire, const AdjacencyMatrix& adjMatrix);

std::vector<double> baseWeights(
    Strategy strategy,
    const Repertoire& repertoire,
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates 
);

std::vector<double> normalizedWeights(
    Strategy strategy,
    const Repertoire& repertoire,
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates
);

// In this version, trait frequency is the probability that a trait is considered for learning
std::vector<double> normalizedWeights(Strategy strategy, const Repertoire& repertoire, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies, std::mt19937& gen);

Repertoire learnTrait(const Repertoire& repertoire, Trait trait);

std::vector<std::pair<Repertoire, double>> transitionFromState(
    Strategy strategy, 
    const Repertoire& repertoire, 
    const AdjacencyMatrix& adjMatrix, 
    const PayoffVector& payoffs, 
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates
);

double stayProbability(std::vector<std::pair<Repertoire, double>> transitions);

struct RepertoireHash {
    std::size_t operator()(const Repertoire& repertoire) const {
        return std::hash<std::string>{}(std::string(repertoire.begin(), repertoire.end()));
    }
};

std::vector<Repertoire> generateReachableRepertoires(
    Strategy strategy, 
    const AdjacencyMatrix& adjMatrix, 
    const PayoffVector& payoffs, 
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates
);

std::vector<Repertoire> generateAllRepertoires(const AdjacencyMatrix& adjMatrix);

size_t countLearnedTraits(const Repertoire& r);

std::vector<Repertoire> retrieveBetterRepertoires(const std::vector<Repertoire>& repertoires, const Repertoire& singleRepertoire);

#endif // LEARNING_HPP