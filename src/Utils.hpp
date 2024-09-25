#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string> 
#include "Types.hpp"

void writeMatrixToCSV(const std::string& filename, const std::vector<std::vector<double>>& matrix);

std::string strategyToString(Strategy strategy);

std::string formatResults(int n, const std::string& adjMatrixBinary, double alpha, Strategy strategy, double expectedSteps);

std::vector<AdjacencyMatrix> readAdjacencyMatrices(int n);

AdjacencyMatrix binaryStringToAdjacencyMatrix(int n, const std::string& str);

bool charToBool(char c);

void printMatrix(const std::vector<std::vector<double>>& matrix);

#endif // UTILS_HPP