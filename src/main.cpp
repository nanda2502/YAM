#include <atomic>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <filesystem>
#include <random>
#include <iomanip>


#include "Debug.hpp"
#include "Utils.hpp"
#include "Types.hpp"
#include "ExpectedSteps.hpp"


void processRepl(
    int repl,
    const AdjacencyMatrix& adjMatrix,
    const Strategy& strategy,
    double alpha,
    const std::vector<size_t>& shuffleSequence,
    int num_steps,
    bool saveTransitionMatrices,
    const std::string &outputDir,
    std::vector<AccumulatedResult>& accumulatedResults,
    size_t idx,
    std::vector<std::atomic<int>>& failureCounts  
) {
    DEBUG_PRINT(1, "Replication:");
    if (DEBUG_LEVEL >= 1) std::cout << repl << '\n';

    std::mt19937 gen(repl);

    double expectedSteps = 0.0;
    double expectedPayoffPerStep = 0.0;
    double expectedTransitionsPerStep = 0.0;
    std::vector<std::vector<double>> transitionMatrix;

    if (!computeExpectedSteps(adjMatrix, strategy, alpha, shuffleSequence, num_steps, expectedSteps, expectedPayoffPerStep, expectedTransitionsPerStep, transitionMatrix)) {
        ++failureCounts[idx];
        return;
    }

    if (saveTransitionMatrices) {
        std::ostringstream alphaStrStream;
        alphaStrStream << std::fixed << std::setprecision(2) << alpha;
        std::string alphaStr = alphaStrStream.str();

        std::string strategyStr = strategyToString(strategy);
        std::string fileName = "transition_mat_" + adjMatrixToBinaryString(adjMatrix) + "_strategy_" + strategyStr + "_alpha_" + alphaStr + ".csv";
        std::string filePath = outputDir + "/" + fileName;
        writeMatrixToCSV(filePath, transitionMatrix);
    }

    // Accumulate results
    AccumulatedResult& accumResult = accumulatedResults[idx];
    accumResult.count++;
    accumResult.totalExpectedSteps += expectedSteps;
    accumResult.totalExpectedTransitionsPerStep += expectedTransitionsPerStep;
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
    try {
        // Define parameters
        std::vector<int> stepVector(20);
        std::iota(stepVector.begin(), stepVector.end(), 1);
        std::vector<double> alphas = {0.0, 0.5};
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
        shuffleSequences.reserve(sequenceCount);  // Reserve space in advance

        do {
            shuffleSequences.push_back(perm);
        } while (std::next_permutation(perm.begin(), perm.end())); 
    
        std::cout << "Starting " << alphas.size() * strategies.size() * adjacencyMatrices.size() * replications * stepVector.size() * shuffleSequences.size() << " runs." << '\n';
    
        // Prepare the combinations
        std::vector<ParamCombination> combinations = makeCombinations(adjacencyMatrices, strategies, alphas, replications, stepVector, shuffleSequences);

        // Accumulate results for each parameter combination excluding shuffle sequence
        std::vector<AccumulatedResult> accumulatedResults(combinations.size());
        std::vector<std::atomic<int>> failureCounts(combinations.size());

        // Parallelize the loop with OpenMP
        #pragma omp parallel for
        for (size_t idx = 0; idx < combinations.size(); ++idx) {
            const ParamCombination& comb = combinations[idx];
            processRepl(
                comb.repl,
                comb.adjMatrix,
                comb.strategy,
                comb.alpha,
                comb.shuffleSequence,
                comb.steps,
                saveTransitionMatrices,
                outputDir,
                accumulatedResults,
                idx,
                failureCounts
            );
        }
    
        // Compute total failures
        int totalFailures = std::accumulate(failureCounts.begin(), failureCounts.end(), 0, [](int sum, const std::atomic<int>& val) {
            return sum + val.load();
        });
    
        // Print total number of failures with debug level 0
        DEBUG_PRINT(0, "Total failures: " << totalFailures);

        // Prepare CSV data with header
        std::string csvHeader = "num_nodes,adj_mat,alpha,strategy,repl,steps,expected_steps,transitions_per_step";
        std::vector<std::string> csvData;
        csvData.push_back(csvHeader);

        // Calculate averages and prepare CSV output
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
                    comb.steps, // Use comb.steps directly since expectedSteps maps to num_steps
                    accumResult.totalExpectedSteps / accumResult.count,
                    accumResult.totalExpectedTransitionsPerStep / accumResult.count
                );
                csvData.push_back(formattedResult);
            }
        }

        // Write and compress CSV output
        writeAndCompressCSV(outputDir, n, csvData);
       
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
    
    return 0;
}