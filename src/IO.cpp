#include "IO.hpp"
#include "Debug.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

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

AdjacencyMatrix parseMatrixString(const std::string& str) {
   // Trim whitespace from the string
   std::string trimmed = str;
   trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);
   
   DEBUG_PRINT(2, "Original string length: " << str.length());
   DEBUG_PRINT(2, "Trimmed string length: " << trimmed.length());
   DEBUG_PRINT(2, "Trimmed string: " << trimmed);
   
   // Check if the trimmed string contains any commas (weighted comma format)
   bool isWeightedCommaFormat = trimmed.find(',') != std::string::npos;
   
   if (isWeightedCommaFormat) {
       // Parse as weighted comma-separated format
       std::vector<double> flatValues;
       std::stringstream ss(trimmed);
       std::string cell;
       
       while (std::getline(ss, cell, ',')) {
           // Handle empty cell (if there are consecutive commas)
           if(cell.empty()) {
               flatValues.push_back(0.0);
           } else {
               try {
                   flatValues.push_back(std::stod(cell));
               } catch (const std::exception& e) {
                   // If conversion fails, default to 0.0
                   flatValues.push_back(0.0);
               }
           }
       }
       
       // Calculate dimension (assuming square matrix)
       int n = static_cast<int>(std::sqrt(flatValues.size()));
       
       // Reshape into n×n matrix
       std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
       for (int i = 0; i < n; i++) {
           for (int j = 0; j < n; j++) {
               size_t index = i * n + j;
               if (index < flatValues.size()) {
                   matrix[i][j] = flatValues[index];
               } else {
                   // Handle case where there are not enough values
                   matrix[i][j] = 0.0;
               }
           }
       }
       
       return matrix;
   } 
   if (trimmed.find_first_not_of("01") == std::string::npos) {
       // Parse as binary format (only contains 0s and 1s)
       int n = static_cast<int>(std::sqrt(trimmed.size()));
       std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
       
       DEBUG_PRINT(2, "Using binary format parser");
       
       for (int i = 0; i < n; i++) {
           for (int j = 0; j < n; j++) {
               // Convert '0'/'1' character to 0.0/1.0 double
               matrix[i][j] = (trimmed[i * n + j] == '1') ? 1.0 : 0.0;
           }
       }
       
       return matrix;
   }
   
   // Parse as compact weighted format (single digits representing tenths)
   DEBUG_PRINT(2, "Using compact weighted format parser");
   int n = static_cast<int>(std::sqrt(trimmed.size()));
   std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
   
   for (int i = 0; i < n; i++) {
       for (int j = 0; j < n; j++) {
           // Convert single digit character to double value with one decimal
           char digit = trimmed[i * n + j];
           int value = digit - '0';  // Convert char to int
           matrix[i][j] = value / 10.0;  // Convert to double with one decimal
       }
   }
   
   return matrix;
}

std::vector<AdjacencyMatrix> readAdjacencyMatrices(const std::string& postfix) {
   std::string filePath = "../data/adj_mat_" + postfix + ".csv";
   std::ifstream file(filePath);
   if (!file.is_open()) {
       throw std::runtime_error("Could not open file " + filePath);
   }

   std::vector<AdjacencyMatrix> matrices;
   std::string line;
   int line_index = 0;
   
   while (std::getline(file, line)) {
       matrices.push_back(parseMatrixString(line));
       
       DEBUG_PRINT(2, "Parsed matrix " << line_index << " from line: " << line);
       if (DEBUG_LEVEL >= 2) {
           std::cout << "Matrix values:" << std::endl;
           const auto& matrix = matrices.back();
           for (const auto & i : matrix) {
               for (double j : i) {
                   std::cout << j << " ";
               }
               std::cout << std::endl;
           }
           std::cout << std::endl;
       }
       
       line_index++;
   }
   
   std::cout << "Loaded " << matrices.size() << " weighted adjacency matrices." << '\n';

   return matrices;
}

void writeAndCompressCSV(const std::string& outputDir, int n, const std::vector<std::string>& csvData) {
    // Construct the output CSV file path (no compression)
    std::string outputCsvPath = outputDir + "/yam_out_" + std::to_string(n) + ".csv";

    std::ofstream csvFile(outputCsvPath);
    if (!csvFile.is_open()) {
        std::cerr << "Failed to open file for writing: " << outputCsvPath << '\n';
        return;
    }
    for (const auto& line : csvData) {
        csvFile << line << "\n";
    }
    csvFile.close();
}

std::string adjMatrixToFlattenedString(const AdjacencyMatrix& adjMatrix) {
    std::stringstream ss;
    size_t n = adjMatrix.size();
    
    // Check if this is a binary matrix (containing only 0.0 and 1.0)
    bool isBinary = true;
    for (size_t i = 0; i < n && isBinary; ++i) {
        for (size_t j = 0; j < n && isBinary; ++j) {
            if (adjMatrix[i][j] != 0.0 && adjMatrix[i][j] != 1.0) {
                isBinary = false;
            }
        }
    }
    
    if (isBinary) {
        // For binary matrices, output as 0s and 1s directly
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ss << (adjMatrix[i][j] == 1.0 ? '1' : '0');
            }
        }
    } else {
        // For weighted matrices, use the compact weighted format
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                // Clamp the value between 0 and 0.9
                double clamped = std::max(0.0, std::min(0.9, adjMatrix[i][j]));
                // Round to 1 decimal place and convert to single digit (0-9)
                int digit = static_cast<int>(std::round(clamped * 10.0));
                if (digit == 10) digit = 9; // Handle potential rounding edge case
                ss << digit;
            }
        }
    }
    
    return ss.str();
}
