#ifndef TYPES_HPP 
#define TYPES_HPP

#include <vector>
#include <string>

using Trait = size_t;
using Repertoire = std::vector<bool>;
using PayoffVector = std::vector<double>;
using AdjacencyMatrix = std::vector<std::vector<bool>>;
using Parents = std::vector<std::vector<Trait>>;

enum Strategy : std::uint8_t {
    Random,
    Payoff,
    Proximal,
    Prestige,
    Conformity,
    Perfect
};

enum traitDistribution : std::uint8_t {
    Learnability,
    Uniform,
    Depth,
    Shallowness,
    Payoffs
};

struct ParamCombination {
    AdjacencyMatrix adjMatrix;
    std::string adjMatrixBinary;
    Strategy strategy;
    traitDistribution distribution;
    double alpha;
    int repl;
    double slope;
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
    double expectedVariation;
};

struct AccumulatedResult {
    int count = 0;
    double absorbing = 0.0;
    std::vector<double> totalExpectedPayoffPerStep{std::vector<double>(20, 0.0)};
    std::vector<double> totalExpectedTransitionsPerStep{std::vector<double>(20, 0.0)};
    std::vector<double> totalExpectedVariation{std::vector<double>(20, 0.0)};
};

#endif // TYPES_HPP

