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
#include <map>

void processRepl(
    int repl,
    const AdjacencyMatrix& adjMatrix,
    const Strategy& strategy,
    double alpha,
    int n,
    double step_factor,
    int num_steps,
    bool saveTransitionMatrices,
    const std::string &outputDir,
    std::vector<Result>& flatResults,
    size_t idx,
    std::vector<std::atomic<int>>& failureCounts  
) {
    DEBUG_PRINT(1, "Replication:");
    if (DEBUG_LEVEL >= 1) std::cout << repl << '\n';

    std::mt19937 gen(repl);

    // Compute expected steps
    double expectedSteps = 0.0;
    double expectedPayoffPerStep = 0.0;
    double expectedTransitionsPerStep = 0.0;
    std::vector<std::vector<double>> transitionMatrix;

    if (!computeExpectedSteps(adjMatrix, strategy, alpha, gen, repl, num_steps, expectedSteps, expectedPayoffPerStep, expectedTransitionsPerStep, transitionMatrix)) {
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
    flatResults[idx] = Result{n, adjMatrixToBinaryString(adjMatrix), alpha, strategy, repl, step_factor, expectedSteps, expectedPayoffPerStep, expectedTransitionsPerStep};
}

std::vector<int> assignSteps(double step_factor, int n, int replications) {
    // Calculate total number of steps
    double step_num = static_cast<double>(n - 1) * step_factor;

    // Determine base steps and fractional steps
    int base_steps = static_cast<int>(std::floor(step_num));
    double frac_steps = step_num - static_cast<double>(base_steps);

    // Calculate how many replications receive an extra step
    int extra_steps_total = static_cast<int>(std::round(frac_steps * static_cast<double>(replications)));

    // Initialize all replications with base_steps
    std::vector<int> num_steps(replications, base_steps);

    if (extra_steps_total > 0) {
        std::vector<int> indices(replications);
        for(int i = 0; i < replications; ++i) {
            indices[i] = i;
        }

        std::random_device rd;
        std::mt19937 gen(rd());

        // Shuffle the indices to randomly select replications for extra steps
        std::shuffle(indices.begin(), indices.end(), gen);

        // Assign an extra step to the first 'extra_steps_total' replications
        for(int i = 0; i < extra_steps_total; ++i) {
            num_steps[indices[i]] += 1;
        }
    }

    return num_steps;
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    bool saveTransitionMatrices = false;
    int numNodes = parseArgs(argc, argv, saveTransitionMatrices);
    int n = numNodes;
    int replications = 30;
    try {
        // Define parameters
        std::vector<double> step_factors = {0.1, 0.25, 0.5, 0.75, 0.9, 1.0, 1.2, 1.5, 2.0};
        std::vector<double> alphas = {1.0};
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
    
        std::cout << "Starting " << alphas.size() * strategies.size() * adjacencyMatrices.size() * replications << " runs." << '\n';
    
        // Prepare the combinations
        std::vector<ParamCombination> combinations = makeCombinations(adjacencyMatrices, strategies, alphas, replications, step_factors);
    
        std::vector<Result> flatResults(combinations.size());
        // Preallocate failure counts vector
        std::vector<std::atomic<int>> failureCounts(combinations.size());
    
        // Create indices for storing the results
        std::vector<size_t> indices(combinations.size());
        std::iota(indices.begin(), indices.end(), 0);
    

        // Precompute steps for each step factor
        std::map<double, std::vector<int>> stepAssignments;
        for (const auto& step_factor : step_factors) {
            stepAssignments[step_factor] = assignSteps(step_factor, n, replications);
        }

        // Process all combinations in parallel
        std::for_each(std::execution::par, indices.begin(), indices.end(), [&](size_t idx) {
            const ParamCombination& comb = combinations[idx];
            DEBUG_PRINT(1, "Adjacency Matrix:");
            if(DEBUG_LEVEL >= 1) std::cout << formatAdjMat(comb.adjMatrixBinary, n) << '\n';
            DEBUG_PRINT(1, "Strategy:");
            if(DEBUG_LEVEL >= 1) std::cout << strategyToString(comb.strategy) << '\n';
            DEBUG_PRINT(1, "Alpha:");
            if(DEBUG_LEVEL >= 1) std::cout << comb.alpha << '\n';

            int num_steps = stepAssignments[comb.step_factor][comb.repl];
    
            processRepl(
                comb.repl,
                comb.adjMatrix,
                comb.strategy,
                comb.alpha,
                n,
                comb.step_factor,
                num_steps,
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
        std::string csvHeader = "num_nodes,adj_mat,alpha,strategy,repl,step_factor,steps,step_payoff,step_transitions";
        std::vector<std::string> csvData;
        csvData.push_back(csvHeader);
    
        for (const auto& result : flatResults) {
            std::string formattedResult = formatResults(
                result.n,
                result.adjMatrixBinary,
                result.alpha,
                result.strategy,
                result.repl,
                result.step_factor,
                result.expectedSteps,
                result.expectedPayoffPerStep,
                result.expectedTransitionsPerStep
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