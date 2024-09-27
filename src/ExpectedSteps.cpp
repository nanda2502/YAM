#include "ExpectedSteps.hpp"

#include "Debug.hpp"
#include "Learning.hpp"
#include "Graph.hpp"
#include "Payoffs.hpp"
#include "LinAlg.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

// Helper function to extract Q matrix and compute (I - Q)
std::vector<std::vector<double>> computeIMinusQ(
    const std::vector<std::vector<double>>& reorderedTransitionMatrix,
    int numTransientStates
) {
    // Extract the Q matrix from the reordered transition matrix. The Q matrix represents the transition probabilities between transient states.
    std::vector<std::vector<double>> qMatrix(numTransientStates, std::vector<double>(numTransientStates));
    for (int i = 0; i < numTransientStates; ++i) {
        for (int j = 0; j < numTransientStates; ++j) {
            qMatrix[i][j] = reorderedTransitionMatrix[i][j];
        }
    }

    DEBUG_PRINT(1, "Q matrix:");
    if(DEBUG_LEVEL >= 1) printMatrix(qMatrix);

    // Subtract the Q matrix from the identity matrix to get I - Q.
    std::vector<std::vector<double>> iMinusQ(numTransientStates, std::vector<double>(numTransientStates));
    for (int i = 0; i < numTransientStates; ++i) {
        for (int j = 0; j < numTransientStates; ++j) {
            iMinusQ[i][j] = (i == j ? 1.0 : 0.0) - qMatrix[i][j];
        }
    }

    DEBUG_PRINT(1, "I - Q matrix:");
    if(DEBUG_LEVEL >= 1) printMatrix(iMinusQ);

    return iMinusQ;
}

std::vector<std::vector<double>> buildTransitionMatrix(
    const std::vector<Repertoire>& repertoiresList,
    const std::unordered_map<Repertoire, int, RepertoireHash>& repertoireIndexMap,
    Strategy strategy,
    const AdjacencyMatrix& adjacencyMatrix,
    const PayoffVector& payoffs
) {
    int numStates = static_cast<int>(repertoiresList.size());
    std::vector<std::vector<double>> transitionMatrix(numStates, std::vector<double>(numStates, 0.0));

    for (int i = 0; i < numStates; ++i) {
        const Repertoire& repertoire = repertoiresList[i];
        auto transitions = transitionFromState(strategy, repertoire, adjacencyMatrix, payoffs);
        double stayProb = stayProbability(strategy, repertoire, adjacencyMatrix, payoffs);

        std::unordered_map<int, double> probMap;
        probMap[i] = stayProb;

        for (const auto& [nextRepertoire, prob] : transitions) {
            int nextIndex = repertoireIndexMap.at(nextRepertoire);
            probMap[nextIndex] += prob;
        }

        std::vector<double> row(numStates, 0.0);
        for (const auto& [idx, prob] : probMap) {
            row[idx] = prob;
        }
        transitionMatrix[i] = row;
    }

    return transitionMatrix;
}

std::tuple<std::vector<std::vector<double>>, std::unordered_map<int, int>, int> reorderTransitionMatrix(
    const std::vector<std::vector<double>>& transitionMatrix,
    const std::vector<std::pair<Repertoire, int>>& repertoiresWithIndices,
    const std::unordered_map<Repertoire, int, RepertoireHash>& repertoireIndexMap,
    Trait rootNode
) {
    // Define the function to check if a state is absorbing
    auto isAbsorbingState = [](const Repertoire& repertoire) {
        return std::all_of(repertoire.begin(), repertoire.end(), [](bool learned) { return learned; });
    };

    // Build a vector of repertoires indexed by state indices
    size_t numStates = transitionMatrix.size();
    std::vector<Repertoire> repertoires(numStates);
    for (const auto& [repertoire, index] : repertoiresWithIndices) {
        repertoires[index] = repertoire;
    }

    // Create a vector of state indices and initialize with 0, 1, ..., numStates - 1
    std::vector<int> reorderedStateIndices(numStates);
    std::iota(reorderedStateIndices.begin(), reorderedStateIndices.end(), 0);

    // Use std::stable_partition to reorder indices: transient states first, then absorbing states
    auto partitionPoint = std::stable_partition(
        reorderedStateIndices.begin(), reorderedStateIndices.end(),
        [&](int index) { return !isAbsorbingState(repertoires[index]); }
    );

    // Calculate the number of transient states
    int numTransientStates = std::distance(reorderedStateIndices.begin(), partitionPoint);

    // Create a mapping from old indices to new indices
    std::unordered_map<int, int> oldToNewIndexMap;
    for (int newIndex = 0; newIndex < static_cast<int>(reorderedStateIndices.size()); ++newIndex) {
        int oldIndex = reorderedStateIndices[newIndex];
        oldToNewIndexMap[oldIndex] = newIndex;
    }

    // Reorder the transition matrix rows and columns
    std::vector<std::vector<double>> reorderedTransitionMatrix;
    reorderedTransitionMatrix.reserve(numStates);

    for (int oldRowIndex : reorderedStateIndices) {
        const auto& oldRow = transitionMatrix[oldRowIndex];
        std::vector<double> newRow;
        newRow.reserve(numStates);
        for (int oldColumnIndex : reorderedStateIndices) {
            newRow.push_back(oldRow[oldColumnIndex]);
        }
        reorderedTransitionMatrix.push_back(std::move(newRow));
    }

    // Identify the initial state index and map it to the new index
    int n = static_cast<int>(repertoiresWithIndices[0].first.size());
    Repertoire initialRepertoire(n, false);
    initialRepertoire[rootNode] = true;

    auto it = repertoireIndexMap.find(initialRepertoire);
    if (it == repertoireIndexMap.end()) {
        throw std::runtime_error("Initial repertoire not found in repertoire index map.");
    }
    int initialStateIndex = it->second;
    int initialStateNewIndex = oldToNewIndexMap.at(initialStateIndex);

    DEBUG_PRINT(1, "Number of transient states: " << numTransientStates);
    DEBUG_PRINT(1, "Number of absorbing states: " << (numStates - numTransientStates));
    DEBUG_PRINT(1, "Initial state new index: " << initialStateNewIndex);

    return {reorderedTransitionMatrix, oldToNewIndexMap, numTransientStates};
}

std::vector<double> computeExpectedPayoffs(
    const std::vector<std::vector<double>>& iMinusQ,
    const std::vector<double>& payoffVector
) {
    // Solve the system (I - Q) * expectedPayoffs = payoffVector
    std::vector<double> expectedPayoffs;
    try {
        expectedPayoffs = solveLinearSystem(iMinusQ, payoffVector);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error computing expected payoffs: " << e.what() << '\n';
        throw;
    }
    return expectedPayoffs;
}

double computeExpectedStepsFromMatrix(
    const std::vector<std::vector<double>>& iMinusQ,
    int initialStateNewIndex
) {
    // Set up the b vector (ones)
    int numTransientStates = static_cast<int>(iMinusQ.size());
    std::vector<double> bVector(numTransientStates, 1.0);

    std::vector<double> tSolution;
    try {
        tSolution = solveLinearSystem(iMinusQ, bVector);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error computing expected steps: " << e.what() << '\n';
        throw;
    }

    double expectedSteps = tSolution[initialStateNewIndex];
    DEBUG_PRINT(1, "Expected steps:");
    if(DEBUG_LEVEL >= 1) std::cout << "Expected steps: " << expectedSteps << '\n';
    return expectedSteps;
}

std::tuple<double, double, std::vector<std::vector<double>>> computeExpectedSteps(
    const AdjacencyMatrix& adjacencyMatrix,
    Strategy strategy,
    double alpha,
    std::mt19937& gen
) {
    Trait rootNode = 0;

    std::vector<int> distances = computeDistances(adjacencyMatrix, rootNode);

    PayoffVector payoffs = generatePayoffs(distances, alpha, gen);

    std::vector<Repertoire> repertoiresList = generateReachableRepertoires(strategy, adjacencyMatrix, payoffs);

    std::vector<std::pair<Repertoire, int>> repertoiresWithIndices;
    for (size_t i = 0; i < repertoiresList.size(); ++i) {
        repertoiresWithIndices.emplace_back(repertoiresList[i], static_cast<int>(i));
    }

    std::unordered_map<Repertoire, int, RepertoireHash> repertoireIndexMap;
    for (const auto& [repertoire, index] : repertoiresWithIndices) {
        repertoireIndexMap[repertoire] = index;
    }

    std::vector<std::vector<double>> transitionMatrix = buildTransitionMatrix(
        repertoiresList, repertoireIndexMap, strategy, adjacencyMatrix, payoffs
    );

    DEBUG_PRINT(1, "Transition matrix:");
    if (DEBUG_LEVEL >= 1) printMatrix(transitionMatrix);

    auto [reorderedTransitionMatrix, oldToNewIndexMap, numTransientStates] = reorderTransitionMatrix(
        transitionMatrix, repertoiresWithIndices, repertoireIndexMap, rootNode
    );

    DEBUG_PRINT(1, "Reordered transition matrix:");
    if (DEBUG_LEVEL >= 1) printMatrix(reorderedTransitionMatrix);

    int initialStateNewIndex = oldToNewIndexMap[repertoireIndexMap[Repertoire(adjacencyMatrix.size(), false)]];
    DEBUG_PRINT(1, "Initial state new index: " << initialStateNewIndex);

    // Compute (I - Q) matrix once
    std::vector<std::vector<double>> iMinusQ = computeIMinusQ(reorderedTransitionMatrix, numTransientStates);

    // Compute expected steps
    double expectedSteps = computeExpectedStepsFromMatrix(
        iMinusQ, initialStateNewIndex
    );

    // Compute state payoffs
    std::vector<double> statePayoffs(repertoiresList.size(), 0.0);
    for (size_t i = 0; i < repertoiresList.size(); ++i) {
        const auto& repertoire = repertoiresList[i];
        for (size_t j = 0; j < repertoire.size(); ++j) {
            if (repertoire[j]) {
                statePayoffs[i] += payoffs[j];
            }
        }
    }

    // Reorder state payoffs according to the new indices
    std::vector<double> reorderedStatePayoffs(statePayoffs.size());
    for (size_t i = 0; i < statePayoffs.size(); ++i) {
        reorderedStatePayoffs[oldToNewIndexMap[i]] = statePayoffs[i];
    }

    // Extract the payoff vector for transient states
    std::vector<double> payoffVector(numTransientStates);
    for (int i = 0; i < numTransientStates; ++i) {
        payoffVector[i] = reorderedStatePayoffs[i];
    }

    // Compute expected payoffs
    std::vector<double> expectedPayoffs = computeExpectedPayoffs(
        iMinusQ, payoffVector
    );

    double totalExpectedPayoff = expectedPayoffs[initialStateNewIndex];
    DEBUG_PRINT(1, "Total expected payoff:");
    if (DEBUG_LEVEL >= 1) std::cout << totalExpectedPayoff << '\n';
    double expectedPayoffPerStep = totalExpectedPayoff / expectedSteps;
    DEBUG_PRINT(1, "Expected payoff per step:");
    if (DEBUG_LEVEL >= 1) std::cout << expectedPayoffPerStep << '\n';

    return {expectedSteps, expectedPayoffPerStep, transitionMatrix};
}

