#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_map>
#include <vector>
#include <string> 
#include "Types.hpp"

void writeMatrixToCSV(const std::string& filename, const std::vector<std::vector<double>>& matrix);

std::string strategyToString(Strategy strategy);

std::string distributionToString(TraitDistribution distribution);

std::string formatResults(
    int n, 
    const std::string& adjMatrixFlattened, 
    double alpha, 
    Strategy strategy, 
    int repl,
    double expectedSteps, 
    double expectedPayoffPerStep, 
    double expectedTransitionsPerStep,
    double expectedVariation,
    double slope,
    TraitDistribution distribution,
    double absorbing,
    double stationaryVariation,
    int payoffDist,
    double edgeWeight,
    double transparency,
    int closure
);

std::vector<AdjacencyMatrix> readAdjacencyMatrices(const std::string& postfix);

std::string formatAdjMat(const std::string& adj_string, int n);

bool charToBool(char c);

void printMatrix(const std::vector<std::vector<double>>& matrix);

int parseArgs(int argc, char* argv[], std::string& postfix);

void writeAndCompressCSV(const std::string& outputDir, int n, const std::vector<std::string>& csvData);

size_t factorial(size_t num);

std::vector<ParamCombination> makeCombinations(
    const std::vector<AdjacencyMatrix>& adjacencyMatrices, 
    int replications,
    const std::string& postfix
);

std::string adjMatrixToFlattenedString(const AdjacencyMatrix& adjMatrix);

std::string stateToString(const Repertoire& state);

void printVector(const std::vector<double>& vec);

void printStates(const std::vector<Repertoire>& repertoiresList, const std::unordered_map<int, int>& oldToNewIndexMap);

AdjacencyMatrix adjustMatrixWeights(const AdjacencyMatrix& adjMatrix, double edgeWeight);

#endif // UTILS_HPP