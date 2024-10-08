#ifndef TYPES_HPP 
#define TYPES_HPP

#include <vector>
#include <string>

using Trait = size_t;
using Repertoire = std::vector<bool>;
using PayoffVector = std::vector<double>;
using AdjacencyMatrix = std::vector<std::vector<bool>>;
using Parents = std::vector<std::vector<Trait>>;

enum Strategy {
    RandomLearning,
    PayoffBasedLearning,
    ProximalLearning,
    PrestigeBasedLearning,
    ConformityBasedLearning
};

struct ParamCombination {
    AdjacencyMatrix adjMatrix;
    std::string adjMatrixBinary;
    Strategy strategy;
    double alpha;
    int repl;
};

struct Result {
    int n;
    std::string adjMatrixBinary;
    double alpha;
    Strategy strategy;
    int repl;
    double expectedSteps;
    double expectedPayoffPerStep;
    double expectedTransitionsPerStep;
};

#endif // TYPES_HPP

