#include <iostream>
#include <mutex>
#include <algorithm>
#include <execution>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <random>
#include <iomanip>
#include <tuple>

#include "Debug.hpp"
#include "Utils.hpp"
#include "Types.hpp"
#include "Graph.hpp"
#include "ExpectedSteps.hpp"

void processRepl(int repl, const AdjacencyMatrix& adjMatrix, const Strategy& strategy, double alpha, int n, bool saveTransitionMatrices, const std::string &outputDir, std::mutex& resultsMutex, std::vector<Result>& flatResults) {
    DEBUG_PRINT(1, "Replication:")
    if (DEBUG_LEVEL >= 1) std::cout << repl << '\n';

    std::random_device rd;
    std::mt19937 gen(rd());

    // Compute expected steps
    double expectedSteps = 0.0;
    double expectedPayoffPerStep = 0.0;
    std::vector<std::vector<double>> transitionMatrix;
    try {
        std::tie(expectedSteps, expectedPayoffPerStep, transitionMatrix) = computeExpectedSteps(adjMatrix, strategy, alpha, gen);
    } catch (const std::exception& e) {
        std::cerr << "Error computing expected steps: " << e.what() << '\n';
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

    // Store results
    std::lock_guard<std::mutex> guard(resultsMutex);
    flatResults.push_back(Result{n, adjMatrixToBinaryString(adjMatrix), alpha, strategy, repl, expectedSteps, expectedPayoffPerStep});
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    bool saveTransitionMatrices = false;
    int numNodes = parseArgs(argc, argv, saveTransitionMatrices);
    int n = numNodes;
    try {
        // Define alphas and strategies
        std::vector<double> alphas = {0.0, 1.0, 2.0};
        std::vector<Strategy> strategies = {Strategy::RandomLearning, Strategy::PayoffBasedLearning};

        // Prepare output directory
        std::string outputDir = "../output";
        if (std::filesystem::exists(outputDir)) {
            std::filesystem::remove_all(outputDir);
        }
        std::filesystem::create_directory(outputDir);

        // Read adjacency matrices
        std::vector<AdjacencyMatrix> adjacencyMatrices = readAdjacencyMatrices(n);

        std::cout << "Loaded " << adjacencyMatrices.size() << " unique adjacency matrices." << '\n';
        std::cout << "Starting " << alphas.size() * strategies.size() * adjacencyMatrices.size() << " runs." << '\n';

        // For each adjacency matrix
        std::vector<Result> flatResults;
        std::mutex resultsMutex;  // Mutex to protect concurrent access to flatResults

        std::for_each(std::execution::par, adjacencyMatrices.begin(), adjacencyMatrices.end(), [&](const auto& adjMatrix) {
            std::string adjMatrixBinary = adjMatrixToBinaryString(adjMatrix);
            DEBUG_PRINT(1, "Adjacency Matrix:");
            if(DEBUG_LEVEL >= 1) std::cout << adjMatrixBinary << '\n';

            std::for_each(std::execution::par, strategies.begin(), strategies.end(), [&](const auto& strategy) {
                DEBUG_PRINT(1, "Strategy:");
                if(DEBUG_LEVEL >= 1) std::cout << strategyToString(strategy) << '\n';

                std::for_each(std::execution::par, alphas.begin(), alphas.end(), [&](const auto& alpha) {
                    DEBUG_PRINT(1, "Alpha:");
                    if(DEBUG_LEVEL >= 1) std::cout << alpha << '\n';

                    std::vector<int> replicas(10);
                    std::iota(replicas.begin(), replicas.end(), 0);

                    std::for_each(std::execution::par, replicas.begin(), replicas.end(), [&](int repl) {
                        processRepl(repl, adjMatrix, strategy, alpha);
                    });
                });
            });
        });

        // Prepare CSV data with header
        std::string csvHeader = "NumNodes,AdjacencyMatrix,Alpha,Strategy,Repl,ExpectedSteps,ExpectedPayoffPerStep";
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

        // Write results to CSV
        std::string outputCsvPath = outputDir + "/expected_steps_" + std::to_string(n) + ".csv";
        std::ofstream csvFile(outputCsvPath);
        if (!csvFile.is_open()) {
            std::cerr << "Failed to open file for writing: " << outputCsvPath << '\n';
            return 1;
        }
        for (const auto& line : csvData) {
            csvFile << line << "\n";
        }
        csvFile.close();

        std::cout << "Expected steps to absorption saved to '" << outputCsvPath << "'" << '\n';

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}