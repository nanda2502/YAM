#include "ExpectedSteps.hpp"

#include "Learning.hpp"
#include "Graph.hpp"
#include "Payoffs.hpp"
#include "LinAlg.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

#define DEBUG_LEVEL 1  // Set to 0 to disable debug output, 1 to enable
#define DEBUG_PRINT(level, x) if (DEBUG_LEVEL >= level) { std::cout << x << '\n'; }

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
    auto isAbsorbingState = [](const Repertoire& repertoire) {
        return std::all_of(repertoire.begin(), repertoire.end(), [](bool learned) { return learned; });
    };

    std::vector<std::pair<int, Repertoire>> absorbingStates;
    std::vector<std::pair<int, Repertoire>> transientStates;

    for (const auto& [repertoire, index] : repertoiresWithIndices) {
        if (isAbsorbingState(repertoire)) {
            absorbingStates.emplace_back(index, repertoire);
        } else {
            transientStates.emplace_back(index, repertoire);
        }
    }

    std::vector<int> reorderedStateIndices;
    reorderedStateIndices.reserve(transientStates.size() + absorbingStates.size());
    for (const auto& [index, _] : transientStates) {
        reorderedStateIndices.push_back(index);
    }
    for (const auto& [index, _] : absorbingStates) {
        reorderedStateIndices.push_back(index);
    }

    std::unordered_map<int, int> oldToNewIndexMap;
    for (int newIndex = 0; newIndex < static_cast<int>(reorderedStateIndices.size()); ++newIndex) {
        int oldIndex = reorderedStateIndices[newIndex];
        oldToNewIndexMap[oldIndex] = newIndex;
    }

    std::vector<std::vector<double>> reorderedTransitionMatrix;
    reorderedTransitionMatrix.reserve(reorderedStateIndices.size());
    for (int oldIndex : reorderedStateIndices) {
        std::vector<double> newRow;
        newRow.reserve(reorderedStateIndices.size());
        for (int oldColumnIndex : reorderedStateIndices) {
            newRow.push_back(transitionMatrix[oldIndex][oldColumnIndex]);
        }
        reorderedTransitionMatrix.push_back(newRow);
    }

    int n = static_cast<int>(repertoiresWithIndices[0].first.size());
    Repertoire initialRepertoire(n, false);
    initialRepertoire[rootNode] = true;

    auto it = repertoireIndexMap.find(initialRepertoire);
    if (it == repertoireIndexMap.end()) {
        throw std::runtime_error("Initial repertoire not found in repertoire index map.");
    }
    int initialStateIndex = it->second;
    int initialStateNewIndex = oldToNewIndexMap[initialStateIndex];

    DEBUG_PRINT(1, "Number of transient states: " << transientStates.size());
    DEBUG_PRINT(1, "Number of absorbing states: " << absorbingStates.size());
    DEBUG_PRINT(1, "Initial state new index: " << initialStateNewIndex);

    return {reorderedTransitionMatrix, oldToNewIndexMap, static_cast<int>(transientStates.size())};
}

double computeExpectedStepsFromMatrix(const std::vector<std::vector<double>>& reorderedTransitionMatrix, int numTransientStates, int initialStateNewIndex) {
    std::cout << "numTransientStates: " << numTransientStates << ", initialStateNewIndex: " << initialStateNewIndex << '\n';
    // Extract the Q matrix from the reordered transition matrix. The Q matrix represents the transition probabilities between transient states.
    std::vector<std::vector<double>> qMatrix(numTransientStates, std::vector<double>(numTransientStates));
    for (int i = 0; i < numTransientStates; ++i) {
        for (int j = 0; j < numTransientStates; ++j) {
            qMatrix[i][j] = reorderedTransitionMatrix[i][j];
        }
    }
    std::cout << "Q matrix:\n";
    printMatrix(qMatrix);

    // Subtract the Q matrix from the identity matrix to get I - Q, representing the needed transformation to find the fundamental matrix.
    std::vector<std::vector<double>> iMinusQ(numTransientStates, std::vector<double>(numTransientStates));
    for (int i = 0; i < numTransientStates; ++i) {
        for (int j = 0; j < numTransientStates; ++j) {
            iMinusQ[i][j] = (i == j ? 1.0 : 0.0) - qMatrix[i][j];
        }
    }

    std::cout << "I - Q matrix:\n";
    printMatrix(iMinusQ);

    std::vector<double> bVector(numTransientStates, 1.0);

    std::vector<double> tSolution;
    try {
        tSolution = solveLinearSystem(iMinusQ, bVector);

    } catch (const std::runtime_error& e) {
        std::cerr << "Error computing expected steps: " << e.what() << '\n';
        throw;
    }

    double expectedSteps = tSolution[initialStateNewIndex];
    std::cout << "Expected steps: " << expectedSteps << '\n';
    return expectedSteps;
}

std::pair<double, std::vector<std::vector<double>>> computeExpectedSteps(const AdjacencyMatrix& adjacencyMatrix, Strategy strategy, double alpha, std::mt19937& gen) {
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

    double expectedSteps = computeExpectedStepsFromMatrix(
        reorderedTransitionMatrix, numTransientStates, initialStateNewIndex
    );

    return {expectedSteps, transitionMatrix};
}

