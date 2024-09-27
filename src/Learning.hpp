#ifndef LEARNING_HPP
#define LEARNING_HPP

#include "Types.hpp"
#include <vector>
#include <string>
#include <random>

std::vector<bool> learnability(const Repertoire& repertoire, const AdjacencyMatrix& adjMatrix);

std::vector<double> baseWeights(Strategy strategy, const PayoffVector& payoffs);

std::vector<double> normalizedWeights(Strategy strategy, const Repertoire& repertoire, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies);

// In this version, trait frequency is the probability that a trait is considered for learning
std::vector<double> normalizedWeights(Strategy strategy, const Repertoire& repertoire, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies, std::mt19937& gen);

Repertoire learnTrait(const Repertoire& repertoire, Trait trait);

std::vector<std::pair<Repertoire, double>> transitionFromState(Strategy strategy, const Repertoire& repertoire, const AdjacencyMatrix& adjMatrix, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies);

double stayProbability(Strategy strategy, const Repertoire& repertoire, const AdjacencyMatrix& adjMatrix, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies);


struct RepertoireHash {
    std::size_t operator()(const Repertoire& repertoire) const {
        return std::hash<std::string>{}(std::string(repertoire.begin(), repertoire.end()));
    }
};

std::vector<Repertoire> generateReachableRepertoires(Strategy strategy, const AdjacencyMatrix& adjMatrix, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies);




#endif // LEARNING_HPP