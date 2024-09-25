#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <filesystem>
#include <random>
#include <stdexcept>
#include <iomanip>
#include <tuple>


#include "Utils.hpp"
#include "Types.hpp"
#include "Graph.hpp"
#include "ExpectedSteps.hpp"

int parseArgs(int argc, char* argv[], bool& saveTransitionMatrices) {
    int numNodes = 4;  // Default value
    saveTransitionMatrices = false;  // Default value

    if (argc == 1) {
        // Use default values
    } else if (argc == 2) {
        numNodes = std::stoi(argv[1]);
        saveTransitionMatrices = true;
    } else if (argc == 3) {
        numNodes = std::stoi(argv[1]);
        std::string arg2 = argv[2];
        if (arg2 == "True") {
            saveTransitionMatrices = true;
        } else if (arg2 == "False") {
            saveTransitionMatrices = false;
        } else {
            throw std::invalid_argument("Second argument must be True or False.");
        }
    } else {
        throw std::invalid_argument("Usage: program [numNodes] [saveTransitionMatrices]");
    }

    return numNodes;
}


int main(int argc, char* argv[])
{
    try {
        // Parse command line arguments
        bool saveTransitionMatrices = false;
        int numNodes = parseArgs(argc, argv, saveTransitionMatrices);
        int n = numNodes;

        // Define alphas and strategies
        std::vector<double> alphas = {0.0, 1.0, 2.0};
        std::vector<Strategy> strategies = {Strategy::RandomLearning, Strategy::PayoffBasedLearning};

        // Initialize random number generator with a fixed seed
        std::mt19937 gen(42);

        // Prepare output directory
        std::string outputDir = "../output";
        if (std::filesystem::exists(outputDir)) {
            std::filesystem::remove_all(outputDir);
        }
        std::filesystem::create_directory(outputDir);

        // Read adjacency matrices
        std::vector<AdjacencyMatrix> adjacencyMatrices = readAdjacencyMatrices(n);

        std::cout << "Loaded " << adjacencyMatrices.size() << " unique adjacency matrices." << '\n';

        // For each adjacency matrix
        std::vector<Result> flatResults;

        for (const auto& adjMatrix : adjacencyMatrices) {
            std::string adjMatrixBinary = adjMatrixToBinaryString(adjMatrix);

            for (const auto& strategy : strategies) {
                for (const auto& alpha : alphas) {
                    // Compute expected steps
                    double expectedSteps = 0.0;
                    std::vector<std::vector<double>> transitionMatrix;
                    try {
                        std::tie(expectedSteps, transitionMatrix) = computeExpectedSteps(adjMatrix, strategy, alpha, gen);
                    } catch (const std::exception& e) {
                        std::cerr << "Error computing expected steps: " << e.what() << '\n';
                        continue;
                    }

                    // Save transition matrix if flag is set
                    if (saveTransitionMatrices) {
                        std::ostringstream alphaStrStream;
                        alphaStrStream << std::fixed << std::setprecision(2) << alpha;
                        std::string alphaStr = alphaStrStream.str();

                        std::string strategyStr = strategyToString(strategy);
                        std::string fileName = "transition_mat_" + adjMatrixBinary + "_strategy_" + strategyStr + "_alpha_" + alphaStr + ".csv";
                        std::string filePath = outputDir + "/" + fileName;
                        writeMatrixToCSV(filePath, transitionMatrix);
                    }

                    // Store results
                    flatResults.push_back(Result{n, adjMatrixBinary, alpha, strategy, expectedSteps});
                }
            }
        }

        // Prepare CSV data with header
        std::string csvHeader = "NumNodes,AdjacencyMatrix,Alpha,Strategy,ExpectedSteps";
        std::vector<std::string> csvData;
        csvData.push_back(csvHeader);

        for (const auto& result : flatResults) {
            std::string formattedResult = formatResults(
                result.n,
                result.adjMatrixBinary,
                result.alpha,
                result.strategy,
                result.expectedSteps
            );
            csvData.push_back(formattedResult);
        }

        // Write results to CSV
        std::string outputCsvPath = outputDir + "/expected_steps.csv";
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