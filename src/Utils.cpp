#include "Utils.hpp"
#include "StringUtils.hpp"
#include <iomanip>
#include <iostream>
#include <stdexcept>

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

int parseArgs(int argc, char* argv[], std::string& postfix) {
    int adj_idx = 1;  // Default value
    postfix = "8";  // Default value

    if (argc == 1) {
        // Use default values
    } else if (argc == 2) {
        adj_idx = std::stoi(argv[1]);
    } else if (argc == 3) {
        adj_idx = std::stoi(argv[1]);
        postfix = std::string(argv[2]);
    } else {
        throw std::invalid_argument("Usage: program [adj_idx] [postfix]");
    }

    return adj_idx;
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
