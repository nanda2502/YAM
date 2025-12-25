#ifndef IO_HPP
#define IO_HPP

#include <string>
#include <vector>
#include "Types.hpp"

// File input/output operations for adjacency matrices and results

// Write a matrix to a CSV file
void writeMatrixToCSV(const std::string& filename, const std::vector<std::vector<double>>& matrix);

// Read adjacency matrices from CSV file
std::vector<AdjacencyMatrix> readAdjacencyMatrices(const std::string& postfix);

// Parse a matrix string from CSV line into AdjacencyMatrix
// Supports three formats:
// 1. Binary flattened (only '0' and '1' characters)
// 2. Compact weighted (single digits 0-9 representing 0.0-0.9)
// 3. Comma-separated weighted (double values separated by commas)
AdjacencyMatrix parseMatrixString(const std::string& str);

// Convert adjacency matrix to flattened string representation
// Binary matrices → string of '0' and '1'
// Weighted matrices → compact format with single digits (0-9)
std::string adjMatrixToFlattenedString(const AdjacencyMatrix& adjMatrix);

// Write CSV data to file (without compression)
void writeAndCompressCSV(const std::string& outputDir, int n, const std::vector<std::string>& csvData);

#endif // IO_HPP
