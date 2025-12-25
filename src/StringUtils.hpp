#ifndef STRINGUTILS_HPP
#define STRINGUTILS_HPP

#include <string>
#include "Types.hpp"

// String formatting and conversion utilities

// Convert Strategy enum to string
std::string strategyToString(Strategy strategy);

// Convert TraitDistribution enum to string
std::string distributionToString(TraitDistribution distribution);

// Format simulation results as CSV line
std::string formatResults(
    int n, 
    const std::string& adjMatrixFlattened, 
    double alpha, 
    Strategy strategy, 
    int repl,
    double expectedSteps, 
    double expectedPayoffPerStep, 
    double expectedTransitionsPerStep,
    double expectedVariation,
    double slope,
    TraitDistribution distribution,
    double absorbing,
    double stationaryVariation,
    int payoffDist,
    double edgeWeight,
    double transparency,
    int closure
);

// Format adjacency matrix string with newlines for display
std::string formatAdjMat(const std::string& adj_string, int n);

// Convert repertoire (binary vector) to string of '0' and '1'
std::string stateToString(const Repertoire& state);

#endif // STRINGUTILS_HPP
