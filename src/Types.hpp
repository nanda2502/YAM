#ifndef TYPES_HPP 
#define TYPES_HPP

#include <vector>
#include <string>

using Trait = size_t;
using Repertoire = std::vector<bool>;
using PayoffVector = std::vector<double>;
using AdjacencyMatrix = std::vector<std::vector<bool>>;

enum Strategy {
    RandomLearning,
    PayoffBasedLearning
};

struct Result {
    int n;
    std::string adjMatrixBinary;
    double alpha;
    Strategy strategy;
    int repl;
    double expectedSteps;
    double expectedPayoffPerStep;
};

#endif // TYPES_HPP

