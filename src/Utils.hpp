#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_map>
#include <vector>
#include <string> 
#include "Types.hpp"

// Basic utility functions for general-purpose operations

// Convert character '0' or '1' to boolean
bool charToBool(char c);

// Print matrix to console
void printMatrix(const std::vector<std::vector<double>>& matrix);

// Print vector to console
void printVector(const std::vector<double>& vec);

// Print repertoire states with reordering
void printStates(const std::vector<Repertoire>& repertoiresList, const std::unordered_map<int, int>& oldToNewIndexMap);

// Parse command-line arguments for adjacency matrix index and postfix
int parseArgs(int argc, char* argv[], std::string& postfix);

#endif // UTILS_HPP