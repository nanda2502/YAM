#include "Utils.hpp"
#include "Types.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <unordered_map>

void writeMatrixToCSV(const std::string& filename, const std::vector<std::vector<double>>& matrix) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }

    for (const auto& row : matrix) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << std::fixed << std::setprecision(4) << row[i];
            if (i < row.size() - 1) file << ',';
        }
        file << '\n';
    }
}

std::string strategyToString(Strategy strategy) {
    switch (strategy) {
        case RandomLearning:
            return "RandomLearning";
        case PayoffBasedLearning:
            return "PayoffBasedLearning";
        case ProximalLearning:
            return "ProximalLearning";
        case PrestigeBasedLearning:
            return "PrestigeBasedLearning";
        case ConformityBasedLearning:
            return "ConformityBasedLearning";
        default:
            throw std::invalid_argument("Unknown strategy");
    }
}

std::string formatResults(
    int n, 
    const std::string& adjMatrixBinary, 
    double alpha, 
    Strategy strategy, 
    int repl,
    double expectedSteps, 
    double expectedPayoffPerStep, 
    double expectedTransitionsPerStep,
    double expectedVariation,
    double slope
) {
    std::ostringstream oss;
    oss << n << ',' << 
    adjMatrixBinary << ',' << 
    alpha << ',' << 
    strategyToString(strategy) << ',' << 
    repl << ',' << 
    std::fixed << std::setprecision(4) << expectedSteps << ',' << 
    expectedPayoffPerStep << ',' << 
    expectedTransitionsPerStep << ',' <<
    expectedVariation << ',' <<
    slope;
    return oss.str();
}

AdjacencyMatrix binaryStringToAdjacencyMatrix(const std::string& str) {
    std::string binaryStr = str;

    int n = std::sqrt(binaryStr.size());

    AdjacencyMatrix matrix(n, std::vector<bool>(n));
    for (int row = 0; row < n; ++row) {
        for (int column = 0; column < n; ++column) {
            matrix[row][column] = charToBool(binaryStr[row * n + column]);
        }
    }

    return matrix;
}

std::vector<AdjacencyMatrix> readAdjacencyMatrices(int n) {
    std::string filePath = "../data/adj_mat_" + std::to_string(n) + ".csv";
    std::ifstream file(filePath);
    if (!file.is_open()) throw std::runtime_error("Could not open file " + filePath);

    std::vector<AdjacencyMatrix> matrices;
    std::string line;
    while (std::getline(file, line)) {
        matrices.push_back(binaryStringToAdjacencyMatrix(line));
    }
    
    std::cout << "Loaded " << matrices.size() << " adjacency matrices." << '\n';

    return matrices;
}

bool charToBool(char c) {
    if (c == '0') return false;
    if (c == '1') return true;
    throw std::invalid_argument("Invalid character in adjacency matrix: Expected '0' or '1', got " + std::string(1, c));
}


void printMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(4) << element << " ";
        }
        std::cout << '\n';
    }
}

int parseArgs(int argc, char* argv[], bool& saveTransitionMatrices) {
    int numNodes = 4;  // Default value
    saveTransitionMatrices = false;  // Default value

    if (argc == 1) {
        // Use default values
    } else if (argc == 2) {
        numNodes = std::stoi(argv[1]);
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

void writeCSV(const std::string& outputDir, int n, const std::vector<std::string>& csvData) {
    // Construct the output CSV file path
    std::string outputCsvPath = outputDir + "/expected_steps_" + std::to_string(n) + ".csv";

    // Write results to CSV
    std::ofstream csvFile(outputCsvPath);
    if (!csvFile.is_open()) {
        std::cerr << "Failed to open file for writing: " << outputCsvPath << '\n';
        return;
    }
    for (const auto& line : csvData) {
        csvFile << line << "\n";
    }
    csvFile.close();

    std::cout << "Expected steps to absorption saved to '" << outputCsvPath << "'\n";
}

std::string adjMatrixToBinaryString(const AdjacencyMatrix& adjMatrix) {
    std::string binaryString;
    binaryString.reserve(adjMatrix.size() * adjMatrix[0].size());

    for (const auto& row : adjMatrix) {
        for (bool entry : row) {
            binaryString += entry ? '1' : '0';
        }
    }
    return binaryString;
}

std::vector<double> returnSlopeVector(Strategy strategy) {
    switch (strategy) {
        case PayoffBasedLearning:
            return {0.0, 1.0,  5.0, 9.0};	
        case ProximalLearning:
            return {1.0,  2.0, 3.0, 5.0};
        case PrestigeBasedLearning:
            return {1.0,  2.0, 3.0,  5.0};
        case ConformityBasedLearning:
            return {0.0, 1.0, 5.0, 15.0};
        default:
            return {0.0};
    }
}

/*
std::vector<double> returnSlopeVector(Strategy strategy) {
    switch (strategy) {
        case PayoffBasedLearning:
            return {5.0};	
        case ProximalLearning:
            return {2.0};
        case PrestigeBasedLearning:
            return {2.0};
        case ConformityBasedLearning:
            return {5.0};
        default:
            return {0.0};
    }
}
*/
std::vector<ParamCombination> makeCombinations(
    const std::vector<AdjacencyMatrix>& adjacencyMatrices, 
    const std::vector<Strategy>& strategies, 
    const std::vector<double>& alphas,
    int replications, 
    const std::vector<int>& stepVector
) {
    std::vector<ParamCombination> combinations;
    for (const auto& adjMatrix : adjacencyMatrices) {
        std::string adjMatrixBinary = adjMatrixToBinaryString(adjMatrix);
        for (const auto& strategy : strategies) {
            auto slopes = returnSlopeVector(strategy);
            for (const auto& alpha : alphas) {
                for (int repl = 0; repl < replications; ++repl) {
                    for (const auto& steps : stepVector) {
                        for (const auto& slope : slopes) {
                        combinations.push_back({adjMatrix, adjMatrixBinary, strategy, alpha, repl, steps, slope});
                        }
                    }
                }
            }
        }
    }
    return combinations;
}

std::string stateToString(const Repertoire& state) {
    std::string binaryString;
    for (bool value : state) {
        binaryString += (value ? '1' : '0');
    }
    return binaryString;
}

void printVector(const std::vector<double>& vec) {
    for (double value : vec) {
        std::cout << value << ' ';
    }
    std::cout << '\n';
}

void printStates(const std::vector<Repertoire>& repertoiresList, const std::unordered_map<int, int>& oldToNewIndexMap) {
    // Generate reordered list of repertoires
    std::vector<Repertoire> reorderedRepertoires(repertoiresList.size());
    for (size_t i = 0; i < repertoiresList.size(); ++i) {
        reorderedRepertoires[oldToNewIndexMap.at(i)] = repertoiresList[i];
    }

    // Print the states (repertoires)
    for (const Repertoire& repertoire : reorderedRepertoires) {
        std::cout << stateToString(repertoire) << '\n';
    }
}