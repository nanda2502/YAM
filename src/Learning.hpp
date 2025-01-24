#ifndef LEARNING_HPP
#define LEARNING_HPP

#include "Types.hpp"
#include <unordered_map>
#include <vector>
#include <string>
#include <random>

std::vector<bool> learnability(const Repertoire& repertoire, const Parents& parents);


// In this version, trait frequency is the probability that a trait is considered for learning
std::vector<double> normalizedWeights(Strategy strategy, const Repertoire& repertoire, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies, std::mt19937& gen, double slope);

Repertoire learnTrait(const Repertoire& repertoire, Trait trait);

double stayProbability(std::vector<std::pair<Repertoire, double>> transitions);

struct RepertoireHash {
    std::size_t operator()(const Repertoire& repertoire) const {
        return std::hash<std::string>{}(std::string(repertoire.begin(), repertoire.end()));
    }
};

std::pair<std::vector<Repertoire>, std::vector<std::vector<std::pair<Repertoire, double>>>>  generateReachableRepertoires(
    Strategy strategy, 
    const AdjacencyMatrix& adjMatrix, 
    const PayoffVector& payoffs, 
    const std::vector<double>& traitFrequencies,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,  
    const std::vector<Repertoire>& allStates,
    const Parents& parents,
    double slope
);

std::vector<Repertoire> generateAllRepertoires(const AdjacencyMatrix& adjMatrix, const Parents& parents);

size_t countLearnedTraits(const Repertoire& r);

std::vector<Repertoire> retrieveBetterRepertoires(const std::vector<Repertoire>& repertoires, const Repertoire& singleRepertoire);

#endif // LEARNING_HPP