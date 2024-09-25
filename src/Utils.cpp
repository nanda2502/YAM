#include "Utils.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <iostream>

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
        default:
            throw std::invalid_argument("Unknown strategy");
    }
}

std::string formatResults(int n, const std::string& adjMatrixBinary, double alpha, Strategy strategy, double expectedSteps) {
    std::ostringstream oss;
    oss << n << ',' << adjMatrixBinary << ',' << alpha << ',' << strategyToString(strategy) << ',' << std::fixed << std::setprecision(2) << expectedSteps;
    return oss.str();
}

std::vector<AdjacencyMatrix> readAdjacencyMatrices(int n) {
    std::string filePath = "../data/adj_mat_" + std::to_string(n) + ".csv";
    std::ifstream file(filePath);
    if (!file.is_open()) throw std::runtime_error("Could not open file " + filePath);

    std::vector<AdjacencyMatrix> matrices;
    std::string line;
    while (std::getline(file, line)) {
        matrices.push_back(binaryStringToAdjacencyMatrix(n, line));
    }

    return matrices;
}

AdjacencyMatrix binaryStringToAdjacencyMatrix(int n, const std::string& str) {
    std::string binaryStr = str;

    if (binaryStr.length() != static_cast<size_t>(n * n)) {
        throw std::invalid_argument("Invalid length: Expected " + std::to_string(n * n) + 
                                    " bits, got " + std::to_string(binaryStr.length()));
    }

    AdjacencyMatrix matrix(n, std::vector<bool>(n));
    for (int row = 0; row < n; ++row) {
        for (int column = 0; column < n; ++column) {
            matrix[row][column] = charToBool(binaryStr[row * n + column]);
        }
    }

    return matrix;
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