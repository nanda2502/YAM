#include "Utils.hpp"
#include "Types.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <unordered_map>
#include <zlib.h>

void writeMatrixToCSV(const std::string& filename, const std::vector<std::vector<double>>& matrix) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }

    for (const auto& row : matrix) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << std::fixed << std::setprecision(2) << row[i];
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
        case SimilarityBasedLearning:
            return "SimilarityBasedLearning";
        case PerfectLearning:
            return "PerfectLearning";
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

    int n = static_cast<int>(std::sqrt(binaryStr.size()));

    if (n == 0) throw std::invalid_argument("Invalid adjmat string: " + str);

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


std::string formatAdjMat(const std::string& adj_string, int n) {
    std::string adj_mat;
    int col_idx = 0;
    for (size_t i = 0; i < adj_string.size(); i++) {
        if (col_idx == n) {
            adj_mat += '\n';
            col_idx = 0;
        }
        adj_mat += adj_string[i];
        col_idx++;
    }
    return adj_mat;
}

bool charToBool(char c) {
    if (c == '0') return false;
    if (c == '1') return true;
    throw std::invalid_argument("Invalid character in adjacency matrix: Expected '0' or '1', got " + std::string(1, c));
}


void printMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(2) << element << " ";
        }
        std::cout << '\n';
    }
}

int parseArgs(int argc, char* argv[], int& num_nodes) {
    int adj_ind = 1;  // Default value
    num_nodes = 8;  // Default value

    if (argc == 1) {
        // Use default values
    } else if (argc == 2) {
        adj_ind = std::stoi(argv[1]);
    } else if (argc == 3) {
        adj_ind = std::stoi(argv[1]);
        num_nodes = std::stoi(argv[2]);
    } else {
        throw std::invalid_argument("Usage: program [adj_ind] [num_nodes]");
    }

    return adj_ind;
}

void writeAndCompressCSV(const std::string& outputDir, int n, const std::vector<std::string>& csvData) {
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

    // Compress the CSV file using gzip
    std::string compressedFilePath = outputCsvPath + ".gz";
    FILE* source = fopen(outputCsvPath.c_str(), "rb");
    gzFile dest = gzopen(compressedFilePath.c_str(), "wb");
    if ((source == nullptr) || (dest == nullptr)) {
        std::cerr << "Failed to open files for compression\n";
        if (source != nullptr) fclose(source);
        if (dest != nullptr) gzclose(dest);
        return;
    }

    char buffer[8192];
    int bytesRead = 0;
    while ((bytesRead = fread(buffer, 1, sizeof(buffer), source)) > 0) {
        gzwrite(dest, buffer, bytesRead);
    }

    fclose(source);
    gzclose(dest);

    // Remove the original uncompressed file
    if (std::remove(outputCsvPath.c_str()) != 0) {
        std::cerr << "Failed to remove original file: " << outputCsvPath << '\n';
    }

    std::cout << "Expected steps to absorption saved and compressed to '" << compressedFilePath << "'\n";
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
/*
std::vector<double> returnSlopeVector(Strategy strategy) {
    if (strategy == RandomLearning || strategy == PerfectLearning) {
        return {0.0};
    } else {
        std::vector<double> slopeVector;
        for (double i = 0.0; i <= 20.0; i += 1.0) {
            slopeVector.push_back(i);
        }
        return slopeVector;
    }
}

*/
std::vector<double> returnSlopeVector(Strategy strategy) {
    switch (strategy) {
        case RandomLearning:
            return {0.0};
        case PerfectLearning:
            return {0.0};
        default:
            return {0.0, 1.0, 1.25, 2.0 , 2.5, 5.0, 10.0, 20.0, 40.0};	

    }
}

/*
std::vector<double> returnSlopeVector(Strategy strategy) {
    switch (strategy) {
        case RandomLearning:
            return {0.0};
        case PerfectLearning:
            return {0.0};
        default:
            return {2.0};	

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