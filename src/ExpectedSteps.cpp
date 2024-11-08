#include "ExpectedSteps.hpp"

#include "Debug.hpp"
#include "Learning.hpp"
#include "Graph.hpp"
#include "Payoffs.hpp"
#include "LinAlg.hpp"
#include "Types.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <numeric>
#include <cmath>


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

    DEBUG_PRINT(2, "Q matrix:");
    if(DEBUG_LEVEL >= 2) printMatrix(qMatrix);

    // Subtract the Q matrix from the identity matrix to get I - Q.
    std::vector<std::vector<double>> iMinusQ(numTransientStates, std::vector<double>(numTransientStates));
    for (int i = 0; i < numTransientStates; ++i) {
        for (int j = 0; j < numTransientStates; ++j) {
            iMinusQ[i][j] = (i == j ? 1.0 : 0.0) - qMatrix[i][j];
        }
    }

    DEBUG_PRINT(2, "I - Q matrix:");
    if(DEBUG_LEVEL >= 2) printMatrix(iMinusQ);

    return iMinusQ;
}

std::vector<std::vector<double>> buildTransitionMatrix(
    const std::vector<Repertoire>& repertoiresList,
    const std::unordered_map<Repertoire, int, RepertoireHash>& repertoireIndexMap,
    Strategy strategy,
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates,
    const Parents& parents
) {
    int numStates = static_cast<int>(repertoiresList.size());
    std::vector<std::vector<double>> transitionMatrix(numStates, std::vector<double>(numStates, 0.0));

    for (int i = 0; i < numStates; ++i) {
        const Repertoire& repertoire = repertoiresList[i];
        auto transitions = transitionFromState(strategy, repertoire, payoffs, traitFrequencies, allStates, parents);
        auto stayProb = stayProbability(transitions);

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

bool isAbsorbingState(const Repertoire& repertoire, const std::vector<double>& traitFrequencies) {
    if (std::all_of(repertoire.begin(), repertoire.end(), [](bool learned) { return learned; })) {
        return true;
    }
    for (Trait i = 1; i < repertoire.size(); ++i) {
        if (!repertoire[i] && std::abs(traitFrequencies[i]) > 1e-5) {
            return false;
        }
    }
    DEBUG_PRINT(1, "No unlearned traits in state " << stateToString(repertoire) << " with non-zero frequency");
    return true;
}

std::tuple<std::vector<std::vector<double>>, std::unordered_map<int, int>, int> reorderTransitionMatrix(
    const std::vector<std::vector<double>>& transitionMatrix,
    const std::vector<std::pair<Repertoire, int>>& repertoiresWithIndices,
    const std::unordered_map<Repertoire, int, RepertoireHash>& repertoireIndexMap,
    Trait rootNode,
    const std::vector<double>& traitFrequencies
) {

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
        [&](int index) { return !isAbsorbingState(repertoires[index], traitFrequencies); }
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

double computeExpectedPayoffOverNSteps(
    const std::vector<std::vector<double>>& transitionMatrix,
    const std::vector<double>& statePayoffs,
    int initialStateIndex,
    double n // Number of steps to simulate
) {
    int numStates = static_cast<int>(transitionMatrix.size());

    // Initialize the state distribution: start in the initial state
    std::vector<double> stateDistribution(numStates, 0.0);
    stateDistribution[initialStateIndex] = 1.0;

    double totalExpectedPayoff = 0.0;

    DEBUG_PRINT(1, "Starting Computation of Expected Payoff Over " << n << " Steps");

    for (int step = 0; step < n; ++step) {
        // Update state distribution for the next step
        std::vector<double> nextStateDistribution(numStates, 0.0);
        for (int i = 0; i < numStates; ++i) {
            for (int j = 0; j < numStates; ++j) {
                nextStateDistribution[j] += stateDistribution[i] * transitionMatrix[i][j];
            }
        }
        stateDistribution = nextStateDistribution;

        // Compute expected payoff at this step
        // Start accumulating payoffs from step 1 and beyond
        double expectedPayoffThisStep = 0.0;
        for (int i = 0; i < numStates; ++i) {
            expectedPayoffThisStep += stateDistribution[i] * statePayoffs[i];
        }
        totalExpectedPayoff += expectedPayoffThisStep;
        DEBUG_PRINT(1, "Expected payoff at Step " << step << " : " << expectedPayoffThisStep);
    }

    DEBUG_PRINT(1, "Total Expected Payoff: " << totalExpectedPayoff);
    
    // Compute expected payoff per step
    double expectedPayoffPerStep = totalExpectedPayoff / n;
    return expectedPayoffPerStep;
}

double computeExpectedStepsFromMatrix(const std::vector<std::vector<double>>& LU, const std::vector<int>& p,  int initialStateNewIndex) {
    // Set up the b vector (ones)
    int numTransientStates = static_cast<int>(LU.size());
    std::vector<double> bVector(numTransientStates, 1.0);

    // Solve the system (I - Q) * t = b
    std::vector<double> tSolution = solveUsingLU(LU, p, bVector);
    
    if (std::any_of(tSolution.begin(), tSolution.end(), [](double val) { return std::isnan(val) || std::isinf(val); })) {
        DEBUG_PRINT(1, "NaN or Inf found in tSolution");
    }

    double expectedSteps = tSolution[initialStateNewIndex];
    if(DEBUG_LEVEL >= 1) std::cout << "Expected steps: " << expectedSteps << '\n';
    return expectedSteps;
}

void adjustTraitFrequencies(std::vector<double>& traitFrequencies, double adjustmentFactor) {
    // Adjust the frequencies slightly to avoid singularity
    size_t n = traitFrequencies.size();
    for (size_t j = 1; j < n; ++j) {  // Exclude trait 0
        traitFrequencies[j] *= adjustmentFactor;  // Reduce frequencies
    }
    // Ensure frequencies are not zero and normalize
    double sum = std::accumulate(traitFrequencies.begin() + 1, traitFrequencies.end(), 0.0);
    if (sum == 0.0) {
        // If all frequencies become zero, reset to equal distribution
        for (size_t j = 1; j < n; ++j) {
            traitFrequencies[j] = 1.0;
        }
        sum = n - 1.0;
    }
    for (size_t j = 1; j < n; ++j) {
        traitFrequencies[j] /= sum;
    }
}

double computeExpectedTransitionsPerStep(
    const std::vector<std::vector<double>>& transitionMatrix,
    int initialStateIndex,
    double n // Number of steps to simulate
) {
    int numStates = static_cast<int>(transitionMatrix.size());

    // Initialize the state distribution: start in the initial state
    std::vector<double> stateDistribution(numStates, 0.0);
    stateDistribution[initialStateIndex] = 1.0;

    double totalExpectedTransitions = 0.0;

    DEBUG_PRINT(1, "Starting Computation of Expected Transitions Over " << n << " Steps");

    for (int step = 0; step < n; ++step) {
        // Compute expected number of transitions at this step
        double expectedTransitionsThisStep = 0.0;

        for (int i = 0; i < numStates; ++i) {
            // Probability of being in state i at the current step
            double probInStateI = stateDistribution[i];

            // Expected number of transitions from state i is (1 - self-transition probability) * probInStateI
            double selfTransitionProb = transitionMatrix[i][i];
            expectedTransitionsThisStep += probInStateI * (1.0 - selfTransitionProb);
        }

        totalExpectedTransitions += expectedTransitionsThisStep;
        DEBUG_PRINT(1, "Expected transitions at Step " << (step + 1) << " : " << expectedTransitionsThisStep);

        // Update state distribution for the next step
        std::vector<double> nextStateDistribution(numStates, 0.0);
        for (int i = 0; i < numStates; ++i) {
            for (int j = 0; j < numStates; ++j) {
                nextStateDistribution[j] += stateDistribution[i] * transitionMatrix[i][j];
            }
        }
        stateDistribution = nextStateDistribution;
    }

    DEBUG_PRINT(1, "Total Expected Transitions: " << totalExpectedTransitions);
    
    // Compute expected transitions per step
    double expectedTransitionsPerStep = totalExpectedTransitions / n;
    return expectedTransitionsPerStep;
}

bool computeExpectedSteps(
    const AdjacencyMatrix& adjacencyMatrix,
    Strategy strategy,
    double alpha,
    const std::vector<size_t>& shuffleSequence,
    int num_steps, 
    double& expectedSteps,                             
    double& expectedPayoffPerStep,
    double& expectedTransitionsPerStep,                     
    std::vector<std::vector<double>>& transitionMatrix 
) {
    try {
        // Initialization
        Strategy baseStrategy = RandomLearning;
        Trait rootNode = 0;
        std::vector<int> distances = computeDistances(adjacencyMatrix, rootNode);
        PayoffVector payoffs = generatePayoffs(distances, alpha, shuffleSequence);

        // print payoffs
        DEBUG_PRINT(2, "Payoffs:");
        if(DEBUG_LEVEL >= 2) {
            for (size_t i = 0; i < payoffs.size(); ++i) {
                std::cout << "Trait " << i << ": " << payoffs[i] << '\n';
            }
        }
        size_t n = adjacencyMatrix.size();
        std::vector<double> traitFrequencies(n, 1.0);
        traitFrequencies[0] = 1.0; // rootNode trait frequency set to 1

        // Compute all parent sets 
        Parents parents(n);
        for (Trait trait = 0; trait < n; ++trait) {
            parents[trait] = parentTraits(adjacencyMatrix, trait);
        }

        // First pass: build initial transition matrix with initial traitFrequencies
        DEBUG_PRINT(1, "Building initial transition matrix with uniform trait frequencies");
        auto allStates = generateAllRepertoires(adjacencyMatrix, parents);

        // Generate repertoires based on initial traitFrequencies
        std::vector<Repertoire> repertoiresList = generateReachableRepertoires(baseStrategy, adjacencyMatrix, payoffs, traitFrequencies, allStates, parents);
        std::vector<std::pair<Repertoire, int>> repertoiresWithIndices;

        for (size_t i = 0; i < repertoiresList.size(); ++i) {
            repertoiresWithIndices.emplace_back(repertoiresList[i], static_cast<int>(i));
        }

        std::unordered_map<Repertoire, int, RepertoireHash> repertoireIndexMap;
        for (const auto& [repertoire, index] : repertoiresWithIndices) {
            repertoireIndexMap[repertoire] = index;
        }

        // Build initial transition matrix
        std::vector<std::vector<double>> preliminaryTransitionMatrix = buildTransitionMatrix(
            repertoiresList, repertoireIndexMap, baseStrategy, payoffs, traitFrequencies, allStates, parents
        );

        auto [reorderedTransitionMatrix, oldToNewIndexMap, numTransientStates] = reorderTransitionMatrix(
            preliminaryTransitionMatrix, repertoiresWithIndices, repertoireIndexMap, rootNode, traitFrequencies
        );

        DEBUG_PRINT(2, "States:");
        if (DEBUG_LEVEL >= 2) printStates(repertoiresList, oldToNewIndexMap);

        DEBUG_PRINT(2, "Initial transition matrix:");        
        if(DEBUG_LEVEL >= 2) printMatrix(reorderedTransitionMatrix);

        if (numTransientStates == 0) {
            expectedSteps = static_cast<double>(num_steps);
            expectedPayoffPerStep = std::accumulate(payoffs.begin(), payoffs.end(), 0.0);
            expectedTransitionsPerStep = 0.0;
            return true;
        }

        // Compute fundamental matrix
        std::vector<std::vector<double>> iMinusQ = computeIMinusQ(reorderedTransitionMatrix, numTransientStates);
        auto [LU, p] = decomposeLU(iMinusQ);
        std::vector<std::vector<double>> fundamentalMatrix(numTransientStates, std::vector<double>(numTransientStates));

        for (int i = 0; i < numTransientStates; ++i) {
            std::vector<double> e_i(numTransientStates, 0.0);
            e_i[i] = 1.0;
            std::vector<double> column = solveUsingLU(LU, p, e_i);

            for (int j = 0; j < numTransientStates; ++j) {
                fundamentalMatrix[j][i] = column[j];
            }
        }

        DEBUG_PRINT(2, "Preliminary Fundamental matrix:");
        if(DEBUG_LEVEL >= 2) printMatrix(fundamentalMatrix);

        double totalTransientTime = 0.0;
        for (const auto& row : fundamentalMatrix) {
            totalTransientTime += std::accumulate(row.begin(), row.end(), 0.0);
        }
        
        // Update trait frequencies based on expected time in states
        for (Trait trait = 1; trait < n; ++trait) {
            double timeTraitKnown = 0.0;
            for (size_t repertoire = 0; repertoire < repertoiresList.size(); ++repertoire) {
                if (repertoiresList[repertoire][trait]) {
                    int newIndex = oldToNewIndexMap[repertoire];
                    if (newIndex < numTransientStates) {
                        timeTraitKnown += fundamentalMatrix[0][newIndex]; 
                    }
                }
            }
            traitFrequencies[trait] = timeTraitKnown / totalTransientTime;
        }

        // Normalize trait frequencies
        double sum = std::accumulate(traitFrequencies.begin() + 1, traitFrequencies.end(), 0.0);
        if (sum == 0.0) {
            // If all frequencies become zero, reset to equal distribution
            for (size_t j = 1; j < n; ++j) {
                traitFrequencies[j] = 1.0;
            }
            sum = n - 1.0;
        }
        for (size_t j = 1; j < n; ++j) {
            traitFrequencies[j] /= sum;
        }

        DEBUG_PRINT(1, "Updated trait frequencies:");
        if(DEBUG_LEVEL >= 1) {
            for (size_t i = 0; i < traitFrequencies.size(); ++i) {
                std::cout << "Trait " << i << ": " << traitFrequencies[i] << '\n';
            }
        }


        // Second pass: rebuild the transition matrix with updated trait frequencies
        DEBUG_PRINT(1, "Building final transition matrix with updated trait frequencies");
        std::vector<Repertoire> finalRepertoiresList = generateReachableRepertoires(
            strategy, adjacencyMatrix, payoffs, traitFrequencies, allStates, parents
        );
        std::unordered_map<Repertoire, int, RepertoireHash> finalRepertoireIndexMap;

        for (size_t i = 0; i < finalRepertoiresList.size(); ++i) {
            finalRepertoireIndexMap[finalRepertoiresList[i]] = static_cast<int>(i);
        }

        // Build final transition matrix without reordering
        transitionMatrix = buildTransitionMatrix(
            finalRepertoiresList, finalRepertoireIndexMap, strategy, payoffs, traitFrequencies, allStates, parents
        );

        if (transitionMatrix.empty()) {
            expectedSteps = static_cast<double>(num_steps);
            expectedPayoffPerStep = std::accumulate(payoffs.begin(), payoffs.end(), 0.0);
            expectedTransitionsPerStep = 0.0;
            return true;
        }

        // Find initial state index
        Repertoire initialRepertoire(n, false);
        initialRepertoire[rootNode] = true;
        int initialStateIndex = finalRepertoireIndexMap[initialRepertoire];

        // Expected steps is the number of steps we simulate
        expectedSteps = static_cast<double>(num_steps);

        // Compute expected payoff per step
        std::vector<double> statePayoffs(finalRepertoiresList.size(), 0.0);
        for (size_t i = 0; i < finalRepertoiresList.size(); ++i) {
            const auto& repertoire = finalRepertoiresList[i];
            for (size_t j = 0; j < repertoire.size(); ++j) {
                if (repertoire[j]) {
                    statePayoffs[i] += payoffs[j];
                }
            }
        }

        expectedPayoffPerStep = computeExpectedPayoffOverNSteps(
            transitionMatrix,
            statePayoffs,
            initialStateIndex,
            num_steps
        );

        // Compute expected transitions per step
        expectedTransitionsPerStep = computeExpectedTransitionsPerStep(
            transitionMatrix,
            initialStateIndex,
            num_steps
        );

        return true;

    } catch (const std::exception& e) {
        // On encountering an exception, return false
        return false;
    }
}