#include <atomic>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <filesystem>
#include "Debug.hpp"
#include "Utils.hpp"
#include "Types.hpp"
#include "ExpectedSteps.hpp"

void processRepl(
    const AdjacencyMatrix& adjMatrix,
    const Strategy& strategy,
    double alpha,
    double slope,
    const std::vector<std::vector<size_t>>& shuffleSequences,
    AccumulatedResult& accumResult,
    std::atomic<int>& failureCount
) {
    double totalExpectedSteps = 0.0;
    double totalExpectedPayoffPerStep = 0.0;
    double totalExpectedTransitionsPerStep = 0.0;
    int successCount = 0;


    for (const auto& shuffleSequence : shuffleSequences) {
        // Compute expected steps
        double expectedSteps = 0.0;
        double expectedPayoffPerStep = 0.0;
        double expectedTransitionsPerStep = 0.0;
        std::vector<std::vector<double>> transitionMatrix;

        if (computeExpectedSteps(adjMatrix, strategy, alpha, shuffleSequence, slope, expectedSteps, expectedPayoffPerStep, expectedTransitionsPerStep, transitionMatrix)) {
            totalExpectedSteps += expectedSteps;
            totalExpectedPayoffPerStep += expectedPayoffPerStep;
            totalExpectedTransitionsPerStep += expectedTransitionsPerStep;
            successCount++;
        } else {
            failureCount++;
        }

    }

    if (successCount > 0) {
        accumResult.count++;
        accumResult.totalExpectedSteps += totalExpectedSteps / successCount;
        accumResult.totalExpectedPayoffPerStep += totalExpectedPayoffPerStep / successCount;
        accumResult.totalExpectedTransitionsPerStep += totalExpectedTransitionsPerStep / successCount;
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
    int n;
    int adj_int = parseArgs(argc, argv, n);
    int replications = 1;

    try {
        // Define alphas and strategies
        std::vector<double> alphas = {0.0};
        std::vector<Strategy> strategies = {
            Strategy::RandomLearning,
            Strategy::PayoffBasedLearning, 
            Strategy::ProximalLearning, 
            Strategy::PrestigeBasedLearning, 
            Strategy::ConformityBasedLearning,
            Strategy::PerfectLearning
        };
    
        // Prepare output directory
        std::string outputDir = "../output";
        if (!std::filesystem::exists(outputDir)) {
            std::filesystem::create_directory(outputDir);
        }
    
        // Read adjacency matrix
        std::vector<AdjacencyMatrix> adjacencyMatricesAll = readAdjacencyMatrices(n);
        std::vector<AdjacencyMatrix> adjacencyMatrices(1, adjacencyMatricesAll[adj_int]);

        std::vector<size_t> perm(n - 1);
        std::iota(perm.begin(), perm.end(), 0);
        size_t sequenceCount = factorial(n - 1);
        std::vector<std::vector<size_t>> shuffleSequences;
        shuffleSequences.reserve(sequenceCount);  

        do {
            shuffleSequences.push_back(perm);
        } while (std::next_permutation(perm.begin(), perm.end())); 
    
        //std::cout << "Starting " << alphas.size() * strategies.size() * adjacencyMatrices.size() * replications * shuffleSequences.size() * 5 << " runs." << '\n';
    
        std::vector<ParamCombination> combinations = makeCombinations(adjacencyMatrices, strategies, alphas, replications);
        std::vector<AccumulatedResult> accumulatedResults(combinations.size());
        std::vector<std::atomic<int>> failureCounts(combinations.size());
    
        std::vector<size_t> indices(combinations.size());
        std::iota(indices.begin(), indices.end(), 0);
        
        #pragma omp parallel for
        for (size_t i = 0; i < indices.size(); ++i) {
            size_t idx = indices[i];
            const auto& comb = combinations[idx];
            processRepl(
                comb.adjMatrix,
                comb.strategy,
                comb.alpha,
                comb.slope,
                shuffleSequences,
                accumulatedResults[idx],
                failureCounts[idx]  
            );
        }
    
        int totalFailures = std::accumulate(failureCounts.begin(), failureCounts.end(), 0, [](int sum, const std::atomic<int>& val) {
            return sum + val.load();
        });
        
        // Print total number of failures with debug level 0
        DEBUG_PRINT(0, "Total failures: " << totalFailures);
    
        // Prepare CSV data with header
        std::string csvHeader = "num_nodes,adj_mat,alpha,strategy,repl,steps,step_payoff,step_transitions,slope";
        std::vector<std::string> csvData;
        csvData.push_back(csvHeader);
    
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
                    accumResult.totalExpectedSteps,
                    accumResult.totalExpectedPayoffPerStep,
                    accumResult.totalExpectedTransitionsPerStep,
                    comb.slope
                );
                csvData.push_back(formattedResult);
            }
        }

    
        writeAndCompressCSV(outputDir, adj_int, csvData);
       
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
    
    return 0;
}