#include <atomic>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <filesystem>
#include <iomanip>

#include "Debug.hpp"
#include "Utils.hpp"
#include "Types.hpp"
#include "ExpectedSteps.hpp"


void processRepl(
    const AdjacencyMatrix& adjMatrix,
    const Strategy& strategy,
    double alpha,
    int steps,
    double slope,
    const std::vector<std::vector<size_t>>& shuffleSequences,
    AccumulatedResult& accumResult,
    std::atomic<int>& failureCount
) {
    double totalExpectedPayoffPerStep = 0.0;
    double totalExpectedTransitionsPerStep = 0.0;
    double totalExpectedVariation = 0.0;
    int successCount = 0;

    for (const auto& shuffleSequence : shuffleSequences) {
        double expectedSteps = 0.0;
        double expectedPayoffPerStep = 0.0;
        double expectedTransitionsPerStep = 0.0;
        double expectedVariation = 0.0;
        std::vector<std::vector<double>> transitionMatrix;


        if (computeExpectedSteps(adjMatrix, strategy, alpha, shuffleSequence, steps, slope, expectedSteps, expectedPayoffPerStep, expectedTransitionsPerStep, expectedVariation, transitionMatrix)) {
            totalExpectedPayoffPerStep += expectedPayoffPerStep;
            totalExpectedTransitionsPerStep += expectedTransitionsPerStep;
            totalExpectedVariation += expectedVariation;
            successCount++;
        } else {
            ++failureCount;
        }
    }

    if (successCount > 0) {
        accumResult.count++;
        accumResult.totalExpectedPayoffPerStep += totalExpectedPayoffPerStep / successCount;
        accumResult.totalExpectedTransitionsPerStep += totalExpectedTransitionsPerStep / successCount;
        accumResult.totalExpectedVariation += totalExpectedVariation / successCount;
    }
}

size_t factorial(size_t num) {
    size_t result = 1;
    for (size_t i = 2; i <= num; ++i) {
        result *= i;
    }
    return result;
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    bool saveTransitionMatrices = false;
    int numNodes = parseArgs(argc, argv, saveTransitionMatrices);
    int n = numNodes;
    int replications = 1;
    std::cout << std::fixed << std::setprecision(4);
    try {
        // Define parameters
        // To do: It would be massively beneficial for performance to instead loop over steps in the computeExpectedSteps function
        std::vector<int> stepVector(20); 
        std::iota(stepVector.begin(), stepVector.end(), 1);
        std::vector<double> alphas = {0.0};
        std::vector<Strategy> strategies = {
            Strategy::RandomLearning,
            Strategy::PayoffBasedLearning,
            Strategy::ProximalLearning,
            Strategy::PrestigeBasedLearning,
            Strategy::ConformityBasedLearning
        };
    
        // Prepare output directory
        std::string outputDir = "../output";
        if (!std::filesystem::exists(outputDir)) {
            std::filesystem::create_directory(outputDir);
        }
    
        // Read adjacency matrices
        std::vector<AdjacencyMatrix> adjacencyMatrices = readAdjacencyMatrices(n);

        // Generate all possible permutations for n - 1 elements
        std::vector<size_t> perm(n - 1);
        std::iota(perm.begin(), perm.end(), 0);
        size_t sequenceCount = factorial(n - 1);
        std::vector<std::vector<size_t>> shuffleSequences;
        shuffleSequences.reserve(sequenceCount);  

        do {
            shuffleSequences.push_back(perm);
        } while (std::next_permutation(perm.begin(), perm.end()));

        
        std::cout << "Starting " << alphas.size() * strategies.size() * adjacencyMatrices.size() * replications * stepVector.size() * shuffleSequences.size() * 5 << " runs." << '\n';
    
        // Prepare the combinations
        std::vector<ParamCombination> combinations = makeCombinations(adjacencyMatrices, strategies, alphas, replications, stepVector);

        // Accumulate results for each parameter combination excluding shuffle sequence
        std::vector<AccumulatedResult> accumulatedResults(combinations.size());
        // Preallocate failure counts vector
        std::vector<std::atomic<int>> failureCounts(combinations.size());
    
        // Create indices for processing
        std::vector<size_t> indices(combinations.size());
        std::iota(indices.begin(), indices.end(), 0);
    
        // Process all combinations in parallel
        #pragma omp parallel for
        for (size_t i = 0; i < indices.size(); ++i) {
            size_t idx = indices[i];
            const ParamCombination& comb = combinations[idx];
            processRepl(
                comb.adjMatrix,
                comb.strategy,
                comb.alpha,
                comb.steps,
                comb.slope,
                shuffleSequences,
                accumulatedResults[idx],
                failureCounts[idx]
            );
        }
    
        // Compute total failures
        int totalFailures = std::accumulate(failureCounts.begin(), failureCounts.end(), 0, [](int sum, const std::atomic<int>& val) {
            return sum + val.load();
        });
    
        // Print total number of failures with debug level 0
        DEBUG_PRINT(0, "Total failures: " << totalFailures);

        // Prepare CSV data with header
        std::string csvHeader = "num_nodes,adj_mat,alpha,strategy,repl,steps,step_payoff,step_transitions,step_variation,slope";
        std::vector<std::string> csvData;
        csvData.push_back(csvHeader);

        // Prepare CSV output
        for (size_t i = 0; i < accumulatedResults.size(); ++i) {
            const AccumulatedResult& accumResult = accumulatedResults[i];
            if (accumResult.count > 0) {
                const ParamCombination& comb = combinations[i];
                std::string formattedResult = formatResults(
                    n,
                    adjMatrixToBinaryString(comb.adjMatrix),
                    comb.alpha,
                    comb.strategy,
                    comb.repl,
                    comb.steps, 
                    accumResult.totalExpectedPayoffPerStep, 
                    accumResult.totalExpectedTransitionsPerStep,
                    accumResult.totalExpectedVariation,
                    comb.slope 
                );
                csvData.push_back(formattedResult);
            }
        }

        // Write output to CSV
        writeCSV(outputDir, n, csvData);
       
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
    
    return 0;
}