#pragma once

#include "Types.hpp"


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
bool computeMarkovChain(
    const AdjacencyMatrix& adjacencyMatrix,
    Strategy strategy, 
    double alpha, // whether payoffs should increase with number of prerequisites
    const std::vector<size_t>& shuffleSequence,
    double slope,  // strength of strategy bias
    double transparency, //knowledge of the learning constraints
    int payoffDist, // 0: equal spacing between 0 and 2 with mean = 1, 1: one high, remaining low values with mean = 1
    TraitDistribution distribution, 
    std::vector<double>& expectedPayoffPerStep,                     
    std::vector<double>& expectedTransitionsPerStep,                
    std::vector<double>& expectedVariation,
    std::vector<std::vector<double>>& transitionMatrix,
    double& timeToAbsorption,
    double& stationaryVariation
);

