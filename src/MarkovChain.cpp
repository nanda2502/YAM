/*
#include "MarkovChain.hpp"

MarkovChain::MarkovChain(
    const ParamCombination& params, 
    const std::vector<size_t>& shuffleSequence,                      
    std::vector<double>& expectedPayoffPerStep,
    std::vector<double>& expectedTransitionsPerStep,
    std::vector<double>& expectedVariation,                       
    std::vector<std::vector<double>>& transitionMatrix,
    double& timeToAbsorption,
    double& stationaryVariation
) : expectedPayoffPerStep(expectedPayoffPerStep),
    expectedTransitionsPerStep(expectedTransitionsPerStep),
    expectedVariation(expectedVariation),
    timeToAbsorption(timeToAbsorption),
    stationaryVariation(stationaryVariation),
    transitionMatrix(transitionMatrix),
    adjMatrix(params.adjMatrix),
    strategy(params.strategy),
    alpha(params.alpha),
    shuffleSequence(shuffleSequence),
    slope(params.slope),
    transparency(params.transparency),
    payoffDist(params.payoffDist),
    distribution(params.distribution),
    baseStrategy(Random),
    rootNode(0), 
    n(adjMatrix.size())
{
    distances = computeDistances(adjMatrix, rootNode);
    payoffs = generatePayoffs(distances, alpha, shuffleSequence, payoffDist);
    DEBUG_PRINT(2, "Payoffs:");
    if(DEBUG_LEVEL >= 2) {
        for (size_t i = 0; i < payoffs.size(); ++i) {
            std::cout << "Trait " << i << ": " << payoffs[i] << '\n';
        }
    }
    DEBUG_PRINT(1, "Slope: " << slope);
    if (DEBUG_LEVEL >= 1) {
        std::cout << "Adjacency Matrix:" << std::endl;
        for (const auto & i : adjMatrix) {
            for (double j : i) {
                std::cout << j << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    std::vector<double> traitFrequencies(n, 1.0);
    traitFrequencies[0] = 1.0; // rootNode trait frequency set to 1

    // Compute all parent sets
    parents.resize(n);
    for (Trait trait = 0; trait < n; ++trait) {
        parents[trait] = parentTraits(adjMatrix, trait);
    }

    // Generate all reachable repertoires
    allStates = generateAllRepertoires(adjMatrix);
}

void MarkovChain::buildInitialTransitionMatrix() {
    DEBUG_PRINT(1, "Building initial transition matrix with uniform trait frequencies");
    
    // Create initial state frequencies (uniform distribution)
    initialStateFrequencies.clear();
    double uniformFrequency = 1.0 / static_cast<double>(allStates.size());
    for (const auto& state : allStates) {
        initialStateFrequencies[state] = uniformFrequency;
    }
    
    // Won't be used, but the function requires it as an argument
    initialStatePayoffs.assign(n, 0.0);

    // Generate repertoires based on initial traitFrequencies using base strategy
    auto [repertoiresList_temp, allTransitions_temp] = computeRepertoiresAndTransitions(
        baseStrategy, slope, initialStatePayoffs, transparency
    );
    repertoiresList = repertoiresList_temp;
    allTransitions = allTransitions_temp;
    
    // Create repertoires with indices
    repertoiresWithIndices.clear();
    repertoiresWithIndices.reserve(repertoiresList.size());
    for (size_t i = 0; i < repertoiresList.size(); ++i) {
        repertoiresWithIndices.emplace_back(repertoiresList[i], static_cast<int>(i));
    }

    // Build repertoire index map
    repertoireIndexMap.clear();
    for (const auto& [repertoire, index] : repertoiresWithIndices) {
        repertoireIndexMap[repertoire] = index;
    }

    // Build preliminary transition matrix
    preliminaryTransitionMatrix = buildTransitionMatrix(
        repertoiresList, repertoireIndexMap, allTransitions
    );

    // Reorder transition matrix
    auto [reorderedMatrix, oldToNew, numTransient] = reorderTransitionMatrix(
        preliminaryTransitionMatrix, repertoiresWithIndices, repertoireIndexMap, rootNode, traitFrequencies
    );
    reorderedTransitionMatrix = reorderedMatrix;
    oldToNewIndexMap = oldToNew;
    numTransientStates = numTransient;

    DEBUG_PRINT(2, "States:");
    if (DEBUG_LEVEL >= 2) printStates(repertoiresList, oldToNewIndexMap);

    DEBUG_PRINT(2, "Initial transition matrix:");        
    if(DEBUG_LEVEL >= 2) printMatrix(reorderedTransitionMatrix);
}

void MarkovChain::computeFundamentalMatrix() {
    // Extract Q matrix from reordered transition matrix
    iMinusQ = computeIMinusQ(reorderedTransitionMatrix, numTransientStates);
    
    // Perform LU decomposition
    auto [LU_temp, p_temp] = decomposeLU(iMinusQ);
    LU = LU_temp;
    p = p_temp;
    
    // Compute fundamental matrix
    fundamentalMatrix.assign(numTransientStates, std::vector<double>(numTransientStates));

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
}

void MarkovChain::calculateStateFrequencies() {
    // Create a mapping between states and their indices in the fundamental matrix
    transientStateIndices.clear();
    transientStates.clear();
    for (size_t i = 0; i < repertoiresList.size(); ++i) {
        int newIndex = oldToNewIndexMap[i];
        if (newIndex < numTransientStates) {
            transientStateIndices[repertoiresList[i]] = newIndex;
            transientStates.push_back(repertoiresList[i]);
        }
    }
    
    // Calculate state frequencies and store them
    stateFrequencies.clear();
    double totalTransientTime = std::accumulate(fundamentalMatrix[0].begin(), fundamentalMatrix[0].end(), 0.0);

    for (const auto& [state, index] : transientStateIndices) {
        double frequency = fundamentalMatrix[0][index] / totalTransientTime;
        stateFrequencies[state] = frequency;
    }

    // Add absorbing state frequency
    absorbingState.assign(n, 1.0);
    stateFrequencies[absorbingState] = 0.05;

    // Normalize frequencies
    double totalStateFreq = 0.0;
    for (const auto& [state, freq] : stateFrequencies) {
        totalStateFreq += freq;
    }
    for (auto& [state, freq] : stateFrequencies) {
        freq /= totalStateFreq;
    }

    DEBUG_PRINT(2, "State Frequencies:");
    if (DEBUG_LEVEL >= 2) {
        for (const auto& [state, freq] : stateFrequencies) {
            std::cout << "State " << stateToString(state) << ": " << freq << '\n';
        }
    }
}

void MarkovChain::updateTraitFrequencies() {
    // Update trait frequencies based on state frequencies
    for (Trait trait = 1; trait < n; ++trait) {
        double timeTraitKnown = 0.0;
        for (const auto& [state, freq] : stateFrequencies) {
            if (state[trait] == 1.0) {
                timeTraitKnown += freq;
            }
        }
        traitFrequencies[trait] = timeTraitKnown;
    }

    // Apply bias adjustments
    traitFrequencies = biasTraitFrequencies(allStates, stateFrequencies, adjMatrix, payoffs, distribution, rootNode);

    DEBUG_PRINT(2, "Adjusted Trait Frequencies:");
    if (DEBUG_LEVEL >= 2) {
        for (Trait trait = 0; trait < n; ++trait) {
            std::cout << "Trait " << trait << ": " << traitFrequencies[trait] << '\n';
        }
    }
    
    // After adjusting trait frequencies, infer compatible state frequencies
    // This is important for strategies like Proximal and Prestige which depend on state frequencies
    if (strategy == Strategy::Proximal || strategy == Strategy::Prestige) {
        DEBUG_PRINT(1, "Inferring state frequencies from target trait frequencies for Proximal/Prestige strategy");
        inferredStateFrequencies = inferStateFrequencies(allStates, traitFrequencies);
    } else {
        // For other strategies, use the original state frequencies
        inferredStateFrequencies = stateFrequencies;
    }
}

void MarkovChain::buildFinalTransitionMatrix() {
    DEBUG_PRINT(1, "Building final transition matrix with updated trait frequencies");
    
    // Compute state payoffs for all states
    allStatesPayoffs.assign(allStates.size(), 0.0);
    for (size_t i = 0; i < allStates.size(); ++i) {
        const auto& state = allStates[i];
        for (size_t j = 0; j < state.size(); ++j) {
            if (state[j] == 1.0) {
                allStatesPayoffs[i] += payoffs[j];
            }
        }
    }

    // Second pass: rebuild the transition matrix with updated trait frequencies
    auto [finalRepertoires, finalTransitions] = computeRepertoiresAndTransitions(
        strategy, slope, allStatesPayoffs, transparency
    );
    finalRepertoiresList = finalRepertoires;
    finalAllTransitions = finalTransitions;

    // Build final repertoire index map
    finalRepertoireIndexMap.clear();
    for (size_t i = 0; i < finalRepertoiresList.size(); ++i) {
        finalRepertoireIndexMap[finalRepertoiresList[i]] = static_cast<int>(i);
    }

    // Build final transition matrix without reordering
    transitionMatrix = buildTransitionMatrix(
        finalRepertoiresList, finalRepertoireIndexMap, finalAllTransitions
    );
}

void MarkovChain::computeResults() {
    // Find initial state index
    initialRepertoire.assign(n, 0.0);
    initialRepertoire[rootNode] = 1.0;
    initialStateIndex = finalRepertoireIndexMap[initialRepertoire];

    // Compute expected payoff per step
    statePayoffs.assign(finalRepertoiresList.size(), 0.0);
    for (size_t i = 0; i < finalRepertoiresList.size(); ++i) {
        const auto& repertoire = finalRepertoiresList[i];
        for (size_t j = 0; j < repertoire.size(); ++j) {
            if (repertoire[j] == 1.0) {
                statePayoffs[i] += payoffs[j];
            }
        }
    }

    timeToAbsorption = computeExpectedTimeToAbsorption(transitionMatrix, initialStateIndex);

    stationaryVariation = computeStationaryVariation(transitionMatrix, finalRepertoiresList);

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
}
*/