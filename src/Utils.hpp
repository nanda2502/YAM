#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <unordered_map>
#include <vector>
#include <string> 
#include "Types.hpp"

void writeMatrixToCSV(const std::string& filename, const std::vector<std::vector<double>>& matrix);

std::string strategyToString(Strategy strategy);

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
);

std::vector<AdjacencyMatrix> readAdjacencyMatrices(int n);

std::string formatAdjMat(const std::string& adj_string, int n);

bool charToBool(char c);

void printMatrix(const std::vector<std::vector<double>>& matrix);

int parseArgs(int argc, char* argv[], int& num_nodes);

void writeAndCompressCSV(const std::string& outputDir, int n, const std::vector<std::string>& csvData);

std::vector<ParamCombination> makeCombinations(
    const std::vector<AdjacencyMatrix>& adjacencyMatrices, 
    const std::vector<Strategy>& strategies, 
    const std::vector<double>& alphas, 
    int replications, 
    const std::vector<int>& stepVector
);

std::string adjMatrixToBinaryString(const AdjacencyMatrix& adjMatrix);

std::string stateToString(const Repertoire& state);

template <typename T>
void printVector(const std::vector<T>& vec) {
    for (const T& value : vec) {
        std::cout << value << ' ';
    }
    std::cout << '\n';
}

void printStates(const std::vector<Repertoire>& repertoiresList, const std::unordered_map<int, int>& oldToNewIndexMap);

#endif // UTILS_HPP