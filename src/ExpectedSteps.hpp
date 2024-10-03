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
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates
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

double computeExpectedTransitionsPerStep(
    const std::vector<std::vector<double>>& fundamentalMatrix,
    const std::vector<std::vector<double>>& reorderedTransitionMatrix,
    int initialStateNewIndex,
    int numTransientStates,
    double expectedSteps
);

// Main function to compute expected steps and expected payoff per step
bool computeExpectedSteps(
    const AdjacencyMatrix& adjacencyMatrix,
    Strategy strategy,
    double alpha,
    std::mt19937& gen,
    double& expectedSteps,                             // Output parameter for expected steps
    double& expectedPayoffPerStep,                     // Output parameter for expected payoff per step
    double& expectedTransitionsPerStep,                // Output parameter for expected transitions per step
    std::vector<std::vector<double>>& transitionMatrix // Output parameter for the transition matrix
);

#endif // EXPECTEDSTEPS_HPP
