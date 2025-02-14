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
#include <map>


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
    std::vector<std::vector<std::pair<Repertoire, double>>>  allTransitions
) {
    int numStates = static_cast<int>(repertoiresList.size());
    std::vector<std::vector<double>> transitionMatrix(numStates, std::vector<double>(numStates, 0.0));

    for (int i = 0; i < numStates; ++i) {
        auto transitions = allTransitions[i];
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

void computeExpectedPayoffAtNSteps(
    const std::vector<std::vector<double>>& transitionMatrix,
    const std::vector<double>& statePayoffs,
    int initialStateIndex,
    std::vector<double>& expectedPayoffPerStep
) {
    int numStates = static_cast<int>(transitionMatrix.size());

    // Initialize the state distribution: start in the initial state
    std::vector<double> stateDistribution(numStates, 0.0);
    stateDistribution[initialStateIndex] = 1.0;

    DEBUG_PRINT(1, "Initial State Distribution:");
    for (int k = 0; k < numStates; ++k) {
        DEBUG_PRINT(1, "State " << k << ": " << stateDistribution[k]);
    }

    for (int step = 0; step < 20; ++step) {
        // Update state distribution for the next step
        std::vector<double> nextStateDistribution(numStates, 0.0);
        for (int i = 0; i < numStates; ++i) {
            for (int j = 0; j < numStates; ++j) {
                nextStateDistribution[j] += stateDistribution[i] * transitionMatrix[i][j];
            }
        }
        stateDistribution = nextStateDistribution; // Update the state distribution

        DEBUG_PRINT(1, "State Distribution at Step " << step + 1 << ":");
        for (int k = 0; k < numStates; ++k) {
            DEBUG_PRINT(1, "State " << k << ": " << stateDistribution[k]);
        }
        double expectedPayoffAtNSteps = 0.0;
        for (int i = 0; i < numStates; ++i) {
            expectedPayoffAtNSteps += stateDistribution[i] * statePayoffs[i];
        }

        DEBUG_PRINT(1, "Expected payoff at Step " << step << " : " << expectedPayoffAtNSteps);
            expectedPayoffPerStep[step] = expectedPayoffAtNSteps;
    }
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

void computeExpectedTransitionsPerStep(
    const std::vector<std::vector<double>>& transitionMatrix,
    int initialStateIndex,
    std::vector<double>& expectedTransitionsPerStep
) {
    int numStates = static_cast<int>(transitionMatrix.size());

    // Initialize the state distribution: start in the initial state
    std::vector<double> stateDistribution(numStates, 0.0);
    stateDistribution[initialStateIndex] = 1.0;

    for (int step = 0; step < 20; ++step) {
        // Compute expected number of transitions at this step
        double expectedTransitionsThisStep = 0.0;

        for (int i = 0; i < numStates; ++i) {
            // Probability of being in state i at the current step
            double probInStateI = stateDistribution[i];

            // Expected number of transitions from state i is (1 - self-transition probability) * probInStateI
            double selfTransitionProb = transitionMatrix[i][i];
            expectedTransitionsThisStep += probInStateI * (1.0 - selfTransitionProb);
        }

        DEBUG_PRINT(1, "Expected transitions at Step " << (step + 1) << " : " << expectedTransitionsThisStep);

        // Store the number of transitions for this step directly in the output vector
        expectedTransitionsPerStep[step] = expectedTransitionsThisStep;

        // Update state distribution for the next step
        std::vector<double> nextStateDistribution(numStates, 0.0);
        for (int i = 0; i < numStates; ++i) {
            for (int j = 0; j < numStates; ++j) {
                nextStateDistribution[j] += stateDistribution[i] * transitionMatrix[i][j];
            }
        }
        stateDistribution = nextStateDistribution;
    }
}

double computeJaccardDistance(const Repertoire& state1, const Repertoire& state2) {
    int intersectionCount = 0;
    int unionCount = 0;

    for (size_t i = 0; i < state1.size(); ++i) {
        if (state1[i] || state2[i]) {
            ++unionCount;
            if (state1[i] && state2[i]) {
                ++intersectionCount;
            }
        }
    }

    return 1.0 - (static_cast<double>(intersectionCount) / static_cast<double>(unionCount));
}

void computeExpectedVariation(const std::vector<std::vector<double>>& transitionMatrix,
                                const std::vector<Repertoire>& repertoires,
                                std::vector<double>& expectedVariation) {
    size_t numStates = transitionMatrix.size();
    std::vector<double> stateProbabilities(numStates, 0.0);
    stateProbabilities[0] = 1.0; 

    // Evolve state probabilities over num_steps
    for (int step = 0; step < 20; ++step) {
        std::vector<double> nextStateProbabilities(numStates, 0.0);

        for (size_t i = 0; i < numStates; ++i) {
            for (size_t j = 0; j < numStates; ++j) {
                nextStateProbabilities[j] += stateProbabilities[i] * transitionMatrix[i][j];
            }
        }

        stateProbabilities = nextStateProbabilities;

        // Calculate total probability mass for active states
        double totalActiveProbability = 0.0;
        for (size_t i = 0; i < numStates; ++i) {
            if (stateProbabilities[i] > 1e-10) {  // Only consider states with non-negligible probability
                totalActiveProbability += stateProbabilities[i];
            }
        }

        // Normalize state probabilities to sum to 1
        std::vector<double> normalizedProbabilities;
        for (size_t i = 0; i < numStates; ++i) {
            if (stateProbabilities[i] > 1e-10) {
                normalizedProbabilities.push_back(stateProbabilities[i] / totalActiveProbability);
            }
        }

        double totalExpectedJaccardDistance = 0.0;
        size_t numComparisons = 0;

        // Calculate average pairwise Jaccard distance between states
        // weighted by their normalized probabilities
        for (size_t i = 0; i < repertoires.size(); ++i) {
            for (size_t j = i + 1; j < repertoires.size(); ++j) {
                if (stateProbabilities[i] > 1e-10 && stateProbabilities[j] > 1e-10) {
                    double jaccardDistance = computeJaccardDistance(repertoires[i], repertoires[j]);
                    double normalizedProb_i = stateProbabilities[i] / totalActiveProbability;
                    double normalizedProb_j = stateProbabilities[j] / totalActiveProbability;
                    totalExpectedJaccardDistance += normalizedProb_i * normalizedProb_j * jaccardDistance;
                    ++numComparisons;
                }
            }
        }

        expectedVariation[step] = (numComparisons == 0) ? 0.0 : totalExpectedJaccardDistance;
    }
}

size_t countLearnedTraits(const Repertoire& repertoire) {
    return std::count(repertoire.begin(), repertoire.end(), true);
}

std::unordered_map<Repertoire, double, RepertoireHash> reverseStateFrequencies(
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies
) {
    std::unordered_map<Repertoire, double, RepertoireHash> reversedFrequencies;
    
    // Group states by number of learned traits
    std::map<size_t, std::vector<std::pair<Repertoire, double>>> statesByTraitCount;
    size_t totalTraits = 0;
    
    // First pass: group states and find total number of traits
    for (const auto& [state, freq] : stateFrequencies) {
        size_t learnedTraits = countLearnedTraits(state);
        statesByTraitCount[learnedTraits].push_back({state, freq});
        totalTraits = std::max(state.size(), totalTraits);
    }
    
    // Second pass: match states and swap frequencies
    for (size_t i = 0; i <= totalTraits; ++i) {
        size_t complementCount = totalTraits - i;
        
        // Skip if either group is empty
        if (statesByTraitCount.find(i) == statesByTraitCount.end() || 
            statesByTraitCount.find(complementCount) == statesByTraitCount.end()) {
            continue;
        }
        
        auto& lowStates = statesByTraitCount[i];
        auto& highStates = statesByTraitCount[complementCount];
        
        // Match states from both groups and swap their frequencies
        size_t pairCount = std::min(lowStates.size(), highStates.size());
        for (size_t j = 0; j < pairCount; ++j) {
            reversedFrequencies[lowStates[j].first] = highStates[j].second;
            reversedFrequencies[highStates[j].first] = lowStates[j].second;
        }
        
        // Handle any remaining unmatched states (keep original frequency)
        for (size_t j = pairCount; j < lowStates.size(); ++j) {
            reversedFrequencies[lowStates[j].first] = lowStates[j].second;
        }
        for (size_t j = pairCount; j < highStates.size(); ++j) {
            reversedFrequencies[highStates[j].first] = highStates[j].second;
        }
    }
    
    return reversedFrequencies;
}

std::vector<double> reverseTraitFrequencies(
    const std::vector<double>& traitFrequencies,
    const AdjacencyMatrix& adjMatrix,
    Trait rootNode
) {
    // Compute distances from root for all traits
    std::vector<int> distances = computeDistances(adjMatrix, rootNode);
    
    // Find maximum distance from root
    int maxDistance = *std::max_element(distances.begin(), distances.end());
    
    // Group traits by their distance from root
    std::map<int, std::vector<std::pair<Trait, double>>> traitsByDistance;
    
    // Skip root node (index 0) when grouping
    for (Trait trait = 1; trait < traitFrequencies.size(); ++trait) {
        traitsByDistance[distances[trait]].push_back({trait, traitFrequencies[trait]});
    }
    
    // Create reversed frequencies vector, initialize with original frequencies
    std::vector<double> reversedFrequencies = traitFrequencies;
    
    // Match and swap frequencies for traits based on complementary distances
    for (int dist = 1; dist <= maxDistance; ++dist) {
        int complementDist = maxDistance - dist + 1;
        
        // Skip if either distance group is empty
        if (traitsByDistance.find(dist) == traitsByDistance.end() || 
            traitsByDistance.find(complementDist) == traitsByDistance.end()) {
            continue;
        }
        
        auto& nearTraits = traitsByDistance[dist];
        auto& farTraits = traitsByDistance[complementDist];
        
        // Match traits from both groups and swap their frequencies
        size_t pairCount = std::min(nearTraits.size(), farTraits.size());
        for (size_t i = 0; i < pairCount; ++i) {
            reversedFrequencies[nearTraits[i].first] = farTraits[i].second;
            reversedFrequencies[farTraits[i].first] = nearTraits[i].second;
        }
        
        // Keep original frequencies for unmatched traits
        for (size_t i = pairCount; i < nearTraits.size(); ++i) {
            reversedFrequencies[nearTraits[i].first] = nearTraits[i].second;
        }
        for (size_t i = pairCount; i < farTraits.size(); ++i) {
            reversedFrequencies[farTraits[i].first] = farTraits[i].second;
        }
    }
    
    // Normalize the frequencies (excluding root node)
    double sum = std::accumulate(
        reversedFrequencies.begin() + 1,  // Skip root node
        reversedFrequencies.end(),
        0.0
    );
    
    // Ensure no zero frequencies and normalize
    for (size_t i = 1; i < reversedFrequencies.size(); ++i) {
        if (reversedFrequencies[i] == 0.0) {
            reversedFrequencies[i] = 1e-5;
            sum += 1e-5;
        }
    }
    
    // Normalize all non-root frequencies
    for (size_t i = 1; i < reversedFrequencies.size(); ++i) {
        reversedFrequencies[i] /= sum;
    }
    
    // Keep root node frequency unchanged at 1.0
    reversedFrequencies[rootNode] = 1.0;
    
    return reversedFrequencies;
}

bool computeExpectedSteps(
    const AdjacencyMatrix& adjacencyMatrix,
    Strategy strategy,
    double alpha,
    const std::vector<size_t>& shuffleSequence,
    double slope,                          
    std::vector<double>& expectedPayoffPerStep,
    std::vector<double>& expectedTransitionsPerStep,
    std::vector<double>& expectedVariation,                       
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

        // this won't be used, but the function requires it as an argument
        std::unordered_map<Repertoire, double, RepertoireHash> initialStateFrequencies;
        double uniformFrequency = 1.0 / static_cast<double>(allStates.size());
        for (const auto& state : allStates) {
            initialStateFrequencies[state] = uniformFrequency;
        }

        // Generate repertoires based on initial traitFrequencies
        auto [repertoiresList, allTransitions] = generateReachableRepertoires(baseStrategy, adjacencyMatrix, payoffs, traitFrequencies, allStates, parents,slope);
        std::vector<std::pair<Repertoire, int>> repertoiresWithIndices;

        repertoiresWithIndices.reserve(repertoiresList.size());
for (size_t i = 0; i < repertoiresList.size(); ++i) {
            repertoiresWithIndices.emplace_back(repertoiresList[i], static_cast<int>(i));
        }

        std::unordered_map<Repertoire, int, RepertoireHash> repertoireIndexMap;
        for (const auto& [repertoire, index] : repertoiresWithIndices) {
            repertoireIndexMap[repertoire] = index;
        }

        // Build initial transition matrix
        std::vector<std::vector<double>> preliminaryTransitionMatrix = buildTransitionMatrix(
            repertoiresList, repertoireIndexMap, allTransitions
        );

        auto [reorderedTransitionMatrix, oldToNewIndexMap, numTransientStates] = reorderTransitionMatrix(
            preliminaryTransitionMatrix, repertoiresWithIndices, repertoireIndexMap, rootNode, traitFrequencies
        );

        DEBUG_PRINT(2, "States:");
        if (DEBUG_LEVEL >= 2) printStates(repertoiresList, oldToNewIndexMap);

        DEBUG_PRINT(2, "Initial transition matrix:");        
        if(DEBUG_LEVEL >= 2) printMatrix(reorderedTransitionMatrix);

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

        // Create a mapping between states and their indices in the fundamental matrix
        std::unordered_map<Repertoire, int, RepertoireHash> transientStateIndices;
        std::vector<Repertoire> transientStates;
        for (size_t i = 0; i < repertoiresList.size(); ++i) {
            int newIndex = oldToNewIndexMap[i];
            if (newIndex < numTransientStates) {
                transientStateIndices[repertoiresList[i]] = newIndex;
                transientStates.push_back(repertoiresList[i]);
            }
        }
        
        // Calculate state frequencies and store them in an unordered_map
        std::unordered_map<Repertoire, double, RepertoireHash> stateFrequencies;
        double totalTransientTime = std::accumulate(fundamentalMatrix[0].begin(), fundamentalMatrix[0].end(), 0.0);

        for (const auto& [state, index] : transientStateIndices) {
            double frequency = fundamentalMatrix[0][index] / totalTransientTime;
            if (frequency == 0.0) {
                frequency = 0.05;  // Frequency of omniscient states
            }
            stateFrequencies[state] = frequency;
        }

        DEBUG_PRINT(2, "State Frequencies:");
        if (DEBUG_LEVEL >= 2) {
            for (const auto& [state, freq] : stateFrequencies) {
                std::cout << "State " << stateToString(state) << ": " << freq << '\n';
            }
        }

        // Update trait frequencies based on state frequencies
        for (Trait trait = 1; trait < n; ++trait) {
            double timeTraitKnown = 0.0;
            for (const auto& [state, freq] : stateFrequencies) {
                if (state[trait]) {
                    timeTraitKnown += freq;
                }
            }
            traitFrequencies[trait] = timeTraitKnown;
        }

        // Normalize trait frequencies
        double sum = std::accumulate(traitFrequencies.begin() + 1, traitFrequencies.end(), 0.0);

        // Ensure no zero frequencies
        for (size_t j = 1; j < n; ++j) {
            if (traitFrequencies[j] == 0.0) {
                traitFrequencies[j] = 0.05;
                sum += 0.05;
            }
        }

        for (size_t j = 1; j < n; ++j) {
            traitFrequencies[j] /= sum;
        }

        // Add trait frequency for absorbing state
        Repertoire absorbingState(n, true);  // All traits learned
        stateFrequencies[absorbingState] = 1e-5;

        //traitFrequencies = reverseTraitFrequencies(traitFrequencies, adjacencyMatrix, 0);
        //stateFrequencies = reverseStateFrequencies(stateFrequencies);

        // Second pass: rebuild the transition matrix with updated trait frequencies
        DEBUG_PRINT(1, "Building final transition matrix with updated trait frequencies");
        auto [finalRepertoiresList, finalAllTransitions] = generateReachableRepertoires(
            strategy, adjacencyMatrix, payoffs, traitFrequencies, allStates, parents, slope
        );
        std::unordered_map<Repertoire, int, RepertoireHash> finalRepertoireIndexMap;

        for (size_t i = 0; i < finalRepertoiresList.size(); ++i) {
            finalRepertoireIndexMap[finalRepertoiresList[i]] = static_cast<int>(i);
        }

        // Build final transition matrix without reordering
        transitionMatrix = buildTransitionMatrix(
            finalRepertoiresList, finalRepertoireIndexMap, finalAllTransitions
        );

        // Find initial state index
        Repertoire initialRepertoire(n, false);
        initialRepertoire[rootNode] = true;
        int initialStateIndex = finalRepertoireIndexMap[initialRepertoire];

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

        computeExpectedPayoffAtNSteps(
            transitionMatrix,
            statePayoffs,
            initialStateIndex,
            expectedPayoffPerStep
        );

        // Compute expected transitions per step
        computeExpectedTransitionsPerStep(
            transitionMatrix,
            initialStateIndex,
            expectedTransitionsPerStep
        );

        // Compute expected variation in traits
        computeExpectedVariation(
            transitionMatrix, 
            finalRepertoiresList, 
            expectedVariation
        );

        return true;

    } catch (const std::exception& e) {
        // On encountering an exception, return false
        return false;
    }
}