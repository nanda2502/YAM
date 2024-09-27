#ifndef EXPECTEDSTEPS_HPP
#define EXPECTEDSTEPS_HPP

#include "Types.hpp"
#include "Learning.hpp"

#include <random>
#include <vector>
#include <unordered_map>
#include <tuple>

// Helper function to compute (I - Q) matrix
std::vector<std::vector<double>> computeIMinusQ(
    const std::vector<std::vector<double>>& reorderedTransitionMatrix,
    int numTransientStates
);

// Function to build the transition matrix
std::vector<std::vector<double>> buildTransitionMatrix(
    const std::vector<Repertoire>& repertoiresList,
    const std::unordered_map<Repertoire, int, RepertoireHash>& repertoireIndexMap,
    Strategy strategy,
    const AdjacencyMatrix& adjacencyMatrix,
    const PayoffVector& payoffs
);

// Function to reorder the transition matrix, separating transient and absorbing states
std::tuple<std::vector<std::vector<double>>, std::unordered_map<int, int>, int> reorderTransitionMatrix(
    const std::vector<std::vector<double>>& transitionMatrix,
    const std::vector<std::pair<Repertoire, int>>& repertoiresWithIndices,
    const std::unordered_map<Repertoire, int, RepertoireHash>& repertoireIndexMap,
    Trait rootNode
);

// Function to compute expected steps from the (I - Q) matrix
double computeExpectedStepsFromMatrix(
    const std::vector<std::vector<double>>& iMinusQ,
    int initialStateNewIndex
);

// Function to compute expected payoffs from the (I - Q) matrix and payoff vector
std::vector<double> computeExpectedPayoffs(
    const std::vector<std::vector<double>>& iMinusQ,
    const std::vector<double>& payoffVector
);

// Main function to compute expected steps and expected payoff per step
std::tuple<double, double, std::vector<std::vector<double>>> computeExpectedSteps(
    const AdjacencyMatrix& adjacencyMatrix,
    Strategy strategy,
    double alpha,
    std::mt19937& gen
);

#endif // EXPECTEDSTEPS_HPP
