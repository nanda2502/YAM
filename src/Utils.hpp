#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_map>
#include <vector>
#include <string> 
#include "Types.hpp"

void writeMatrixToCSV(const std::string& filename, const std::vector<std::vector<double>>& matrix);

std::string strategyToString(Strategy strategy);

std::string formatResults(int n, const std::string& adjMatrixBinary, double alpha, Strategy strategy, int repl, double expectedSteps, double expectedPayoffPerStep, double expectedTransitionsPerStep);

std::vector<AdjacencyMatrix> readAdjacencyMatrices(int n);

AdjacencyMatrix binaryStringToAdjacencyMatrix(int n, const std::string& str);

std::string formatAdjMat(const std::string& adj_string, int n);

bool charToBool(char c);

void printMatrix(const std::vector<std::vector<double>>& matrix);

int parseArgs(int argc, char* argv[], bool& saveTransitionMatrices);

void writeAndCompressCSV(const std::string& outputDir, int n, const std::vector<std::string>& csvData);

std::string adjMatrixToBinaryString(const AdjacencyMatrix& adjMatrix);

std::vector<ParamCombination> makeCombinations(std::vector<AdjacencyMatrix>& adjacencyMatrices, std::vector<Strategy>& strategies, std::vector<double>& alphas, int replications);

std::string stateToString(const Repertoire& state);

void printStates(const std::vector<Repertoire>& repertoiresList, const std::unordered_map<int, int>& oldToNewIndexMap);

#endif // UTILS_HPP