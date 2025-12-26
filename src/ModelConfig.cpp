#include "ModelConfig.hpp"
#include "MarkovChain.hpp"
#include "IO.hpp"
#include "StringUtils.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <unordered_set>

std::vector<double> returnSlopeVector(Strategy strategy) {
    switch (strategy) {
        case Random: 
            return {0.0};
        default:
            return {1.25, 2.0, 5.0};	
    }
}

bool isUnconstrained(const AdjacencyMatrix& adjMatrix) {
    // Skip the first row
    for (size_t row = 1; row < adjMatrix.size(); ++row) {
        // Check each element in this row
        for (double col : adjMatrix[row]) {
            // If any element is true (nonzero), return false
            if (col > 0.0) {
                return false;
            }
        }
    }
    // If we've checked all rows and found no nonzero values, return true
    return true;
}

size_t factorial(size_t num) {
    size_t result = 1;
    for (size_t i = 2; i <= num; ++i) {
        result *= i;
    }
    return result;
}

std::vector<std::vector<size_t>> makeShuffles(int n) {
    std::vector<std::vector<size_t>> shuffleSequences;
    
    // For small n (10 or less), generate all permutations
    if (n <= 10) {
        std::vector<size_t> perm(n - 1);
        std::iota(perm.begin(), perm.end(), 0);
        size_t sequenceCount = factorial(n - 1);
        shuffleSequences.reserve(sequenceCount);
        
        do {
            shuffleSequences.push_back(perm);
        } while (std::ranges::next_permutation(perm).found);
    } 
    // For larger n, generate a limited number of random permutations
    else {
        // Generate random unique permutations of trait indices (excluding root)
        size_t maxPermutations = 10000;
        
        // Create a set to store unique permutations (as strings for easy comparison)
        std::unordered_set<std::string> uniquePermsSet;
        
        // Create a random number generator
        std::random_device rd;
        std::mt19937 g(rd());
        
        // Generate the base permutation
        std::vector<size_t> basePerm(n - 1);
        std::iota(basePerm.begin(), basePerm.end(), 0);
        
        // Try to generate maxPermutations unique permutations
        while (uniquePermsSet.size() < maxPermutations) {
            // Create a new permutation by shuffling the base
            std::vector<size_t> perm = basePerm;
            std::shuffle(perm.begin(), perm.end(), g);
            
            // Convert to string for uniqueness check
            std::string permStr;
            for (size_t val : perm) {
                permStr += std::to_string(val) + ",";
            }
            
            // Add to set and vector if unique
            if (uniquePermsSet.insert(permStr).second) {
                shuffleSequences.push_back(perm);
                
                // Progress output
                if (shuffleSequences.size() % 500 == 0) {
                    std::cout << "Generated " << shuffleSequences.size() << " unique permutations..." << '\n';
                }
            }
            
            // Safety check to avoid infinite loop
            if (uniquePermsSet.size() == factorial(n - 1)) {
                std::cout << "Generated all possible permutations (" << uniquePermsSet.size() << ")" << '\n';
                break;
            }
        }
    }
    
    return shuffleSequences;
}

AdjacencyMatrix computeTransitiveClosure(const AdjacencyMatrix& adjacencyMatrix) {
    size_t n = adjacencyMatrix.size();
    
    // Make a copy of the original matrix to work with
    AdjacencyMatrix result = adjacencyMatrix;
    
    // Floyd-Warshall algorithm for transitive closure
    for (size_t k = 0; k < n; ++k) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                // If there's a path from i to k and from k to j, then there's a path from i to j
                if (result[i][k] > 0.0 && result[k][j] > 0.0) {
                    result[i][j] = 1.0; // Set to 1.0 to indicate a path exists
                }
            }
        }
    }
    
    return result;
}

// Read this function to understand the specific parameter combinations used for each figure in the paper. 
std::vector<ParamCombination> makeCombinations(
    const std::vector<AdjacencyMatrix>& adjacencyMatrices, 
    int replications,
    const std::string& postfix
) {
    
    constexpr bool SENSITIVITY_TESTS = true; // Set to false to skip sensitivity tests
    std::vector<ParamCombination> combinations;
    
    // Define default values
    TraitDistribution defaultDistribution = TraitDistribution::Learnability;
    int defaultPayoffDist = 0;
    double defaultAlpha = 0.0;
    double alternativeAlpha = 1.0;
    double defaultEdgeWeight = 1.0;
    double defaultTransparency = 0.0;
    int defaultClosure = 0; // By default, dont use transitive closure
    
    std::vector<Strategy> strategies = {
        Strategy::Random,
        Strategy::Payoff,
        Strategy::Proximal,
        Strategy::Prestige,
        Strategy::Conformity,
        Strategy::Prestige2 // Weight by repertoire size instead of repertoire payoff
    };

    std::vector<TraitDistribution> distributions = {
        TraitDistribution::Learnability,
        TraitDistribution::Uniform,
        TraitDistribution::Depth,
        TraitDistribution::Shallowness,
        TraitDistribution::Payoffs
    };

    std::vector<double> weights(21, 0.0);
    for (int i = 0; i < 21; ++i) weights[i] = i * 0.05;

    for (const auto & adjMatrix : adjacencyMatrices) {
         std::string adjMatrixFlattened = adjMatrixToFlattenedString(adjMatrix);
        size_t n = adjMatrix.size();
        auto shuffleSequences = makeShuffles(n);
        // Determine which shuffle sequences to use
        std::vector<std::vector<size_t>> usedShuffleSequences;
        if (isUnconstrained(adjMatrix) || n == 121) {
            // For unconstrained adjacency matrix, use only the first shuffle sequence
            usedShuffleSequences = {shuffleSequences[0]};
            std::cout << "Using only the first shuffle sequence" << '\n';
        } else {
            // For constrained adjacency matrix, use all shuffle sequences
            usedShuffleSequences = shuffleSequences;
        }
        
        // Create combinations with default parameters
        for (const auto& strategy : strategies) {
            // Default slope based on strategy 
            double defaultSlope = strategy == Strategy::Random ? 0.0 : 2.0;

            if (n <= 8) {
                
                // For n <= 8, use all slopes
                auto slopes = returnSlopeVector(strategy);
                // These parameters were used for Figure 1, Figure 2 A, Figure S2-S4, Figure S6, Figure S10. 
                // Slope = 2 was used for everything except Figure S4, which shows slope = 1.25 and 5.0.
                
                // Base cases: default values for all parameters, but vary the slopes
                for (const auto& slope : slopes) {
                    for (int repl = 0; repl < replications; ++repl) {
                        combinations.push_back({
                            adjMatrix, 
                            adjMatrixFlattened, 
                            strategy, 
                            defaultDistribution, 
                            defaultAlpha, 
                            repl, 
                            slope, 
                            defaultPayoffDist, 
                            usedShuffleSequences,
                            defaultEdgeWeight,
                            defaultTransparency,
                            defaultClosure
                        });
                    }
                }
                
                // Only continue with parameter variation if the adjacency matrix is size 8
                if (n == 8 && SENSITIVITY_TESTS) {
                    
                    // Figure 2F
                    // Vary alpha: Alpha controls whether payoffs are randomly distributed (alpha = 0.0) or increase with number of prerequisites (alpha = 1.0)
                    for (int repl = 0; repl < replications; ++repl) {
                        combinations.push_back({
                            adjMatrix, 
                            adjMatrixFlattened, 
                            strategy, 
                            defaultDistribution, 
                            alternativeAlpha, 
                            repl, 
                            defaultSlope, 
                            defaultPayoffDist, 
                            usedShuffleSequences,
                            defaultEdgeWeight,
                            defaultTransparency,
                            defaultClosure
                        });
                    }
                    
                    // Figure 2 D & E
                    // Vary trait expression distribution: advancedness and high payoff expression bias 
                    for (const auto& distribution : distributions) {
                        if (distribution != defaultDistribution) {
                            for (int repl = 0; repl < replications; ++repl) {
                                combinations.push_back({
                                    adjMatrix, 
                                    adjMatrixFlattened, 
                                    strategy, 
                                    distribution, 
                                    defaultAlpha, 
                                    repl, 
                                    defaultSlope, 
                                    defaultPayoffDist, 
                                    usedShuffleSequences,
                                    defaultEdgeWeight,
                                    defaultTransparency,
                                    defaultClosure
                                });
                            }
                        }
                    }
                    
                    // Figure S8
                    // Vary payoffDist: Add combinations with one payoff being high and the rest low
                    for (size_t payoffDist = 1; payoffDist < 2; ++payoffDist) {  // Start from 1 since 0 is default
                        for (int repl = 0; repl < replications; ++repl) {
                            combinations.push_back({
                                adjMatrix, 
                                adjMatrixFlattened, 
                                strategy, 
                                defaultDistribution, 
                                defaultAlpha, 
                                repl, 
                                defaultSlope, 
                                static_cast<int>(payoffDist), 
                                usedShuffleSequences,
                                defaultEdgeWeight,
                                defaultTransparency,
                                defaultClosure
                            });
                        }
                    }
                    // Figure 2 G-I                     
                    // Vary edge weight: we handle the possibility of skipping by converting the structure to its transitive closure
                    // i.e., all ancestors of a trait are considered prerequisites
                    auto transitiveClosure = computeTransitiveClosure(adjMatrix);
                    auto transitiveClosureFlattened = adjMatrixToFlattenedString(transitiveClosure);

                    for (const auto& weight : weights){
                        for (int repl = 0; repl < replications; ++repl) {
                            combinations.push_back({
                                transitiveClosure, 
                                transitiveClosureFlattened, 
                                strategy, 
                                defaultDistribution, 
                                defaultAlpha, 
                                repl, 
                                defaultSlope, 
                                defaultPayoffDist, 
                                usedShuffleSequences,
                                weight,
                                defaultTransparency,
                                1
                            });
                        }
                    }
                    
                    // Figure S12                     
                    // Vary edge weight: transitive reduction
                    // i.e., no ancestors of a trait are considered prerequisites

                    for (const auto& weight : weights){
                        for (int repl = 0; repl < replications; ++repl) {
                            combinations.push_back({
                                adjMatrix, 
                                adjMatrixFlattened, 
                                strategy, 
                                defaultDistribution, 
                                defaultAlpha, 
                                repl, 
                                defaultSlope, 
                                defaultPayoffDist, 
                                usedShuffleSequences,
                                weight,
                                defaultTransparency,
                                defaultClosure
                            });
                        }
                    }
                    
                    
                    // Figure 2 J-L
                    // Add combinations with different transparency values
                    for (double transparency : {0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0}) {
                        for (int repl = 0; repl < replications; ++repl) {
                            combinations.push_back({
                                adjMatrix, 
                                adjMatrixFlattened, 
                                strategy, 
                                defaultDistribution, 
                                defaultAlpha, 
                                repl, 
                                defaultSlope, 
                                defaultPayoffDist, 
                                usedShuffleSequences,
                                defaultEdgeWeight,
                                transparency,
                                defaultClosure
                            });
                        }
                    }
                    
                    
                }

            } else {
                // For n > 8, use only the default slope
                double defaultSlope = strategy == Strategy::Random ? 0.0 : 2.0;
                
                for (int repl = 0; repl < replications; ++repl) {
                    combinations.push_back({
                        adjMatrix, 
                        adjMatrixFlattened, 
                        strategy, 
                        defaultDistribution, 
                        defaultAlpha, 
                        repl, 
                        defaultSlope, 
                        defaultPayoffDist, 
                        usedShuffleSequences,
                        defaultEdgeWeight,
                        defaultTransparency,
                        defaultClosure
                    });
                }
            }
        }
    }
    
    return combinations;
}

AdjacencyMatrix adjustMatrixWeights(const AdjacencyMatrix& adjMatrix, double edgeWeight) {
    auto weightedMatrix = adjMatrix;
    for (size_t row = 0; row < adjMatrix.size(); row++) {
        for (size_t col = 0; col < adjMatrix.size(); col++) {
            // if there's an edge, adjust the edge weight
            weightedMatrix[row][col] = adjMatrix[row][col] == 1.0 ? edgeWeight : 0.0;
        }
    }
    return weightedMatrix;
}

void processRepl(
    const ParamCombination& params,
    AccumulatedResult& accumResult,
    std::atomic<int>& failureCount
) {
    std::cout << "Processing strategy: " << strategyToString(params.strategy) << '\n';
    std::vector<double> totalExpectedPayoffPerStep(20, 0.0);
    std::vector<double> totalExpectedTransitionsPerStep(20, 0.0);
    double totalTimeToAbsorption = 0.0;
    int successCount = 0;

    AdjacencyMatrix weightedMatrix;
    if (params.edgeWeight != 1.0) {
        weightedMatrix = adjustMatrixWeights(params.adjMatrix, params.edgeWeight);
    } else {
        weightedMatrix = params.adjMatrix;
    }

    for (const auto& shuffleSequence : params.shuffleSequences) {
        double timeToAbsorption;
        timeToAbsorption = shuffleSequence == params.shuffleSequences[0] ? 0.0 : -1.0; // only compute time to absorption for the first shuffle sequence
        std::vector<double> expectedPayoffPerStep(20, 0.0);
        std::vector<double> expectedTransitionsPerStep(20, 0.0);
        std::vector<std::vector<double>> transitionMatrix;

        if (computeMarkovChain(weightedMatrix, params.strategy, params.alpha, shuffleSequence,
                                 params.slope, params.transparency, params.payoffDist, params.distribution, 
                                 expectedPayoffPerStep,
                                 expectedTransitionsPerStep,
                                 transitionMatrix, 
                                 timeToAbsorption)) {
            totalTimeToAbsorption += timeToAbsorption;
            for (size_t i = 0; i < 20; ++i) {
                totalExpectedPayoffPerStep[i] += expectedPayoffPerStep[i];
                totalExpectedTransitionsPerStep[i] += expectedTransitionsPerStep[i];
            }
          successCount++;
        } else {
          failureCount++;
        }
    }

    if (successCount > 0) {
        for (size_t i = 0; i < 20; ++i) {
            accumResult.count++;
            accumResult.totalExpectedPayoffPerStep[i] += totalExpectedPayoffPerStep[i] / params.shuffleSequences.size();
            accumResult.totalExpectedTransitionsPerStep[i] += totalExpectedTransitionsPerStep[i] / params.shuffleSequences.size();
        }
        accumResult.absorbing += totalTimeToAbsorption / params.shuffleSequences.size();
    }

    std::cout << "Failures: " << failureCount.load() << '\n';
}
