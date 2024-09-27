#include <iostream>
#include <algorithm>
#include <execution>
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
    int n,
    bool saveTransitionMatrices,
    const std::string &outputDir,
    std::vector<Result>& flatResults,
    size_t idx,
    std::vector<std::atomic<int>>& failureCounts  // Add failureCounts vector
) {
    DEBUG_PRINT(1, "Replication:");
    if (DEBUG_LEVEL >= 1) std::cout << repl << '\n';

    std::random_device rd;
    std::mt19937 gen(rd());

    // Compute expected steps
    double expectedSteps = 0.0;
    double expectedPayoffPerStep = 0.0;
    std::vector<std::vector<double>> transitionMatrix;

    if (!computeExpectedSteps(adjMatrix, strategy, alpha, gen, expectedSteps, expectedPayoffPerStep, transitionMatrix)) {
        // Increment failure count
        ++failureCounts[idx];
        return;
    }

    // Save transition matrix if flag is set
    if (saveTransitionMatrices) {
        std::ostringstream alphaStrStream;
        alphaStrStream << std::fixed << std::setprecision(2) << alpha;
        std::string alphaStr = alphaStrStream.str();

        std::string strategyStr = strategyToString(strategy);
        std::string fileName = "transition_mat_" + adjMatrixToBinaryString(adjMatrix) + "_strategy_" + strategyStr + "_alpha_" + alphaStr + ".csv";
        std::string filePath = outputDir + "/" + fileName;
        writeMatrixToCSV(filePath, transitionMatrix);
    }

    // Store results directly into the pre-allocated vector
    flatResults[idx] = Result{n, adjMatrixToBinaryString(adjMatrix), alpha, strategy, repl, expectedSteps, expectedPayoffPerStep};
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    bool saveTransitionMatrices = false;
    int numNodes = parseArgs(argc, argv, saveTransitionMatrices);
    int n = numNodes;
    int replications = 10;
    try {
        // Define alphas and strategies
        std::vector<double> alphas = {0.0, 1.0, 2.0};
        std::vector<Strategy> strategies = {Strategy::RandomLearning, Strategy::PayoffBasedLearning};
    
        // Prepare output directory
        std::string outputDir = "../output";
        if (!std::filesystem::exists(outputDir)) {
            std::filesystem::create_directory(outputDir);
        }
    
        // Read adjacency matrices
        std::vector<AdjacencyMatrix> adjacencyMatrices = readAdjacencyMatrices(n);
    
        std::cout << "Starting " << alphas.size() * strategies.size() * adjacencyMatrices.size() * 10 << " runs." << '\n';
    
        // Prepare the combinations
        std::vector<ParamCombination> combinations = makeCombinations(adjacencyMatrices, strategies, alphas, replications);
    
        std::vector<Result> flatResults(combinations.size());
        // Preallocate failure counts vector
        std::vector<std::atomic<int>> failureCounts(combinations.size());
    
        // Create indices for storing the results
        std::vector<size_t> indices(combinations.size());
        std::iota(indices.begin(), indices.end(), 0);
    
        // Process all combinations in parallel
        std::for_each(std::execution::par, indices.begin(), indices.end(), [&](size_t idx) {
            const ParamCombination& comb = combinations[idx];
            DEBUG_PRINT(1, "Adjacency Matrix:");
            if(DEBUG_LEVEL >= 1) std::cout << comb.adjMatrixBinary << '\n';
            DEBUG_PRINT(1, "Strategy:");
            if(DEBUG_LEVEL >= 1) std::cout << strategyToString(comb.strategy) << '\n';
            DEBUG_PRINT(1, "Alpha:");
            if(DEBUG_LEVEL >= 1) std::cout << comb.alpha << '\n';
    
            processRepl(
                comb.repl,
                comb.adjMatrix,
                comb.strategy,
                comb.alpha,
                n,
                saveTransitionMatrices,
                outputDir,
                flatResults,
                idx,
                failureCounts  
            );
        });
    
        // Compute total failures
        int totalFailures = 0;
        for (const auto& count : failureCounts) {
            totalFailures += count.load();
        }
    
        // Print total number of failures with debug level 0
        DEBUG_PRINT(0, "Total failures: " << totalFailures);
    
        // Prepare CSV data with header
        std::string csvHeader = "num_nodes,adj_mat,alpha,strategy,repl,steps,step_payoff";
        std::vector<std::string> csvData;
        csvData.push_back(csvHeader);
    
        for (const auto& result : flatResults) {
            std::string formattedResult = formatResults(
                result.n,
                result.adjMatrixBinary,
                result.alpha,
                result.strategy,
                result.repl,
                result.expectedSteps,
                result.expectedPayoffPerStep
            );
            csvData.push_back(formattedResult);
        }
    
        writeAndCompressCSV(outputDir, n, csvData);
       
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
    
    return 0;
}