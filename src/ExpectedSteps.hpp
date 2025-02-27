#ifndef EXPECTEDSTEPS_HPP
#define EXPECTEDSTEPS_HPP

#include "Types.hpp"
#include "Learning.hpp"

#include <vector>
#include <unordered_map>
#include <tuple>


std::vector<std::vector<double>> computeIMinusQ(
    const std::vector<std::vector<double>>& reorderedTransitionMatrix,
    int numTransientStates
);

std::vector<std::vector<double>> buildTransitionMatrix(
    const std::vector<Repertoire>& repertoiresList,
    const std::unordered_map<Repertoire, int, RepertoireHash>& repertoireIndexMap,
    std::vector<std::vector<std::pair<Repertoire, double>>>  allTransitions
);

bool isAbsorbingState(const Repertoire& repertoire, const std::vector<double>& traitFrequencies);

// Function to reorder the transition matrix, separating transient and absorbing states
std::tuple<std::vector<std::vector<double>>, std::unordered_map<int, int>, int> reorderTransitionMatrix(
    const std::vector<std::vector<double>>& transitionMatrix,
    const std::vector<std::pair<Repertoire, int>>& repertoiresWithIndices,
    const std::unordered_map<Repertoire, int, RepertoireHash>& repertoireIndexMap,
    Trait rootNode,
    const std::vector<double>& traitFrequencies
);


double computeExpectedStepsFromMatrix(
    const std::vector<std::vector<double>>& LU,
    const std::vector<int>& p,
    int initialStateNewIndex
);

// Solve the linear system using the payoff vector
std::vector<double> computeExpectedPayoffs(
    const std::vector<std::vector<double>>& LU,
    const std::vector<int>& p,
    const std::vector<double>& payoffVector
);

// Compute learning success probability
double computeExpectedTransitionsPerStep(
    const std::vector<std::vector<double>>& fundamentalMatrix,
    const std::vector<std::vector<double>>& reorderedTransitionMatrix,
    int initialStateNewIndex,
    int numTransientStates,
    double expectedSteps
);

// Main function
bool computeExpectedSteps(
    const AdjacencyMatrix& adjacencyMatrix,
    Strategy strategy,
    double alpha,
    const std::vector<size_t>& shuffleSequence,
    int num_steps,
    double slope, 
    double& expectedSteps,                             // Output parameter for expected steps
    double& expectedPayoffPerStep,                     // Output parameter for expected payoff per step
    double& expectedTransitionsPerStep,                // Output parameter for expected transitions per step
    double& expectedVariation,                          // Output parameter for expected variation in traits
    std::vector<std::vector<double>>& transitionMatrix // Output parameter for the transition matrix
);

#endif // EXPECTEDSTEPS_HPP
