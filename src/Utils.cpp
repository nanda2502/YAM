#include "Utils.hpp"
#include "Types.hpp"
#include "Debug.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
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
        case Random:
            return "Random";
        case Payoff:
            return "Payoff";
        case Proximal:
            return "Proximal";
        case Prestige:
            return "Prestige";
        case Conformity:
            return "Conformity";
        case Perfect:
            return "Perfect";
        case Anticonformity:
            return "Anticonformity";
        case Prestige2:
            return "Prestige2";
        default:
            throw std::invalid_argument("Unknown strategy");
    }
}

std::string distributionToString(TraitDistribution distribution) {
    switch (distribution) {
        case Learnability:
            return "Learnability";
        case Uniform:
            return "Uniform";
        case Depth:
            return "Depth";
        case Shallowness:
            return "Shallowness";
        case Payoffs:
            return "Payoffs";
        default:
            throw std::invalid_argument("Unknown distribution");
    }
}

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
) {
    std::ostringstream oss;
    oss << n << ',' << 
    adjMatrixFlattened << ',' << 
    alpha << ',' << 
    strategyToString(strategy) << ',' << 
    repl << ',' << 
    std::fixed << std::setprecision(4) << expectedSteps << ',' << 
    expectedPayoffPerStep << ',' << 
    expectedTransitionsPerStep << ',' <<
    expectedVariation << ',' <<
    slope << ',' <<
    distributionToString(distribution) << ',' <<
    absorbing << ',' <<
    stationaryVariation << ',' <<
    payoffDist << ',' <<
    edgeWeight << ',' <<
    transparency << ',' <<
    closure;
    return oss.str();
}

inline std::vector<std::vector<double>> flattenedStringToWeightedMatrix(const std::string& str) {
    // For backward compatibility, convert flattened string to weighted matrix
    // where 1s become 1.0 and 0s become 0.0
    
    int n = static_cast<int>(std::sqrt(str.size()));
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
    
    for (int row = 0; row < n; ++row) {
        for (int column = 0; column < n; ++column) {
            // Convert the character '0' or '1' to double 0.0 or 1.0
            matrix[row][column] = (str[row * n + column] == '1') ? 1.0 : 0.0;
        }
    }
    
    return matrix;
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

std::vector<double> returnSlopeVector(Strategy strategy) {
    switch (strategy) {
        case Random: case Perfect:
            return {0.0};
        default:
            return {1.25, 2.0, 5.0};	

    }
}

/*
std::vector<double> returnSlopeVector(Strategy strategy) {
    switch (strategy) {
        case Random:
            return {0.0};
        case Perfect:
            return {0.0};
        default:
            return {2.0};	

    }
}
*/
bool isUnconstrained(const AdjacencyMatrix& adjMatrix) {
    // Skip the first row
    for (size_t row = 1; row < adjMatrix.size(); ++row) {
        // Check each element in this row
        for (double col : adjMatrix[row]) {
            // If any element is true (nonzero), return false
            if (col > 0.0) {
                return false;
            }
        }
    }
    // If we've checked all rows and found no nonzero values, return true
    return true;
}


size_t factorial(size_t num) {
    size_t result = 1;
    for (size_t i = 2; i <= num; ++i) {
        result *= i;
    }
    return result;
}

std::vector<std::vector<size_t>> makeShuffles(int n) {
    std::vector<std::vector<size_t>> shuffleSequences;
    
    // For small n (10 or less), generate all permutations
    if (n <= 10) {
        std::vector<size_t> perm(n - 1);
        std::iota(perm.begin(), perm.end(), 0);
        size_t sequenceCount = factorial(n - 1);
        shuffleSequences.reserve(sequenceCount);
        
        do {
            shuffleSequences.push_back(perm);
        } while (std::ranges::next_permutation(perm).found);
    } 
    // For larger n, generate a limited number of random permutations
    else {
        // Generate random unique permutations of trait indices (excluding root)
        size_t maxPermutations = 10000;
        
        // Create a set to store unique permutations (as strings for easy comparison)
        std::unordered_set<std::string> uniquePermsSet;
        
        // Create a random number generator
        std::random_device rd;
        std::mt19937 g(rd());
        
        // Generate the base permutation
        std::vector<size_t> basePerm(n - 1);
        std::iota(basePerm.begin(), basePerm.end(), 0);
        
        // Try to generate maxPermutations unique permutations
        while (uniquePermsSet.size() < maxPermutations) {
            // Create a new permutation by shuffling the base
            std::vector<size_t> perm = basePerm;
            std::shuffle(perm.begin(), perm.end(), g);
            
            // Convert to string for uniqueness check
            std::string permStr;
            for (size_t val : perm) {
                permStr += std::to_string(val) + ",";
            }
            
            // Add to set and vector if unique
            if (uniquePermsSet.insert(permStr).second) {
                shuffleSequences.push_back(perm);
                
                // Progress output
                if (shuffleSequences.size() % 500 == 0) {
                    std::cout << "Generated " << shuffleSequences.size() << " unique permutations..." << '\n';
                }
            }
            
            // Safety check to avoid infinite loop
            if (uniquePermsSet.size() == factorial(n - 1)) {
                std::cout << "Generated all possible permutations (" << uniquePermsSet.size() << ")" << '\n';
                break;
            }
        }
    }
    
    return shuffleSequences;
}

AdjacencyMatrix computeTransitiveClosure(const AdjacencyMatrix& adjacencyMatrix) {
    size_t n = adjacencyMatrix.size();
    
    // Make a copy of the original matrix to work with
    AdjacencyMatrix result = adjacencyMatrix;
    
    // Floyd-Warshall algorithm for transitive closure
    for (size_t k = 0; k < n; ++k) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                // If there's a path from i to k and from k to j, then there's a path from i to j
                if (result[i][k] > 0.0 && result[k][j] > 0.0) {
                    result[i][j] = 1.0; // Set to 1.0 to indicate a path exists
                }
            }
        }
    }
    
    return result;
}

// Read this function to understand the specific parameter combinations used for each figure in the paper. 
std::vector<ParamCombination> makeCombinations(
    const std::vector<AdjacencyMatrix>& adjacencyMatrices, 
    int replications,
    const std::string& postfix
) {
    constexpr bool SENSITIVITY_TESTS = false; // Set to false to skip sensitivity tests
    std::vector<ParamCombination> combinations;
    
    // Define default values
    TraitDistribution defaultDistribution = TraitDistribution::Learnability;
    int defaultPayoffDist = 0;
    double defaultAlpha = 0.0;
    double alternativeAlpha = 1.0;
    double defaultEdgeWeight = 1.0;
    double defaultTransparency = 0.0;
    int defaultClosure = 0; // By default, dont use transitive closure
    
    std::vector<Strategy> strategies = {
        //Strategy::Random,
        //Strategy::Payoff,
        Strategy::Proximal//,
        //Strategy::Prestige,
        //Strategy::Conformity,
        //Strategy::Perfect,
        //Strategy::Anticonformity,
        //Strategy::Prestige2
    };

    std::vector<TraitDistribution> distributions = {
        TraitDistribution::Learnability,
        TraitDistribution::Uniform,
        TraitDistribution::Depth,
        TraitDistribution::Shallowness,
        TraitDistribution::Payoffs
    };

    std::vector<double> weights(21, 0.0);
    for (int i = 0; i < 21; ++i) weights[i] = i * 0.05;

    for (const auto & adjMatrix : adjacencyMatrices) {
         std::string adjMatrixFlattened = adjMatrixToFlattenedString(adjMatrix);
        size_t n = adjMatrix.size();
        auto shuffleSequences = makeShuffles(n);
        // Determine which shuffle sequences to use
        std::vector<std::vector<size_t>> usedShuffleSequences;
        if (isUnconstrained(adjMatrix) || n == 121) {
            // For unconstrained adjacency matrix, use only the first shuffle sequence
            usedShuffleSequences = {shuffleSequences[0]};
            std::cout << "Using only the first shuffle sequence" << '\n';
        } else {
            // For constrained adjacency matrix, use all shuffle sequences
            usedShuffleSequences = shuffleSequences;
        }
        
        // Create combinations with default parameters
        for (const auto& strategy : strategies) {
            // Default slope based on strategy 
            double defaultSlope = (strategy == Strategy::Random || strategy == Strategy::Perfect) ? 0.0 : 2.0;

            if (n <= 8) {
                
                // For n <= 8, use all slopes
                auto slopes = returnSlopeVector(strategy);
                // These parameters were used for figure 1, figure 2 E-H, A, Q. Slope 2 was used for everything except figure 2 E-H, which shows slope 1.25 and 5.0.
                // Base cases: default values for all parameters, but vary the slopes
                for (const auto& slope : slopes) {
                    for (int repl = 0; repl < replications; ++repl) {
                        combinations.push_back({
                            adjMatrix, 
                            adjMatrixFlattened, 
                            strategy, 
                            defaultDistribution, 
                            defaultAlpha, 
                            repl, 
                            slope, 
                            defaultPayoffDist, 
                            usedShuffleSequences,
                            defaultEdgeWeight,
                            defaultTransparency,
                            defaultClosure
                        });
                    }
                }
                
                // Only continue with parameter variation if the adjacency matrix is size 8
                if (n == 8 && SENSITIVITY_TESTS) {
                    /*
                    // These parameters were used for figure 2 L
                    // Vary alpha: Alpha controls whether payoffs are randomly distributed (alpha = 0.0) or increase with number of prerequisites (alpha = 1.0)
                    for (int repl = 0; repl < replications; ++repl) {
                        combinations.push_back({
                            adjMatrix, 
                            adjMatrixFlattened, 
                            strategy, 
                            defaultDistribution, 
                            alternativeAlpha, 
                            repl, 
                            defaultSlope, 
                            defaultPayoffDist, 
                            usedShuffleSequences,
                            defaultEdgeWeight,
                            defaultTransparency,
                            defaultClosure
                        });
                    }
                    
                    // Figure 2 D & E
                    // Vary distribution: Add combinations with each non-default distribution
                    for (const auto& distribution : distributions) {
                        if (distribution != defaultDistribution) {
                            for (int repl = 0; repl < replications; ++repl) {
                                combinations.push_back({
                                    adjMatrix, 
                                    adjMatrixFlattened, 
                                    strategy, 
                                    distribution, 
                                    defaultAlpha, 
                                    repl, 
                                    defaultSlope, 
                                    defaultPayoffDist, 
                                    usedShuffleSequences,
                                    defaultEdgeWeight,
                                    defaultTransparency,
                                    defaultClosure
                                });
                            }
                        }
                    }
                    
                    // Figure SX
                    // Vary payoffDist: Add combinations with one payoff being high and the rest low
                    for (size_t payoffDist = 1; payoffDist < 2; ++payoffDist) {  // Start from 1 since 0 is default
                        for (int repl = 0; repl < replications; ++repl) {
                            combinations.push_back({
                                adjMatrix, 
                                adjMatrixFlattened, 
                                strategy, 
                                defaultDistribution, 
                                defaultAlpha, 
                                repl, 
                                defaultSlope, 
                                static_cast<int>(payoffDist), 
                                usedShuffleSequences,
                                defaultEdgeWeight,
                                defaultTransparency,
                                defaultClosure
                            });
                        }
                    }
                    // Figure 2 G-I                     
                    // Vary edge weight: we handle the possibility of skipping by converting the structure to its transitive closure
                    // i.e., all ancestors of a trait are considered prerequisites
                    auto transitiveClosure = computeTransitiveClosure(adjMatrix);
                    auto transitiveClosureFlattened = adjMatrixToFlattenedString(transitiveClosure);

                    for (const auto& weight : weights){
                        for (int repl = 0; repl < replications; ++repl) {
                            combinations.push_back({
                                transitiveClosure, 
                                transitiveClosureFlattened, 
                                strategy, 
                                defaultDistribution, 
                                defaultAlpha, 
                                repl, 
                                defaultSlope, 
                                defaultPayoffDist, 
                                usedShuffleSequences,
                                weight,
                                defaultTransparency,
                                1
                            });
                        }
                    }
                    
                    // Figure SX                     
                    // Vary edge weight: transitive reduction
                    // i.e., no ancestors of a trait are considered prerequisites

                    for (const auto& weight : weights){
                        for (int repl = 0; repl < replications; ++repl) {
                            combinations.push_back({
                                adjMatrix, 
                                adjMatrixFlattened, 
                                strategy, 
                                defaultDistribution, 
                                defaultAlpha, 
                                repl, 
                                defaultSlope, 
                                defaultPayoffDist, 
                                usedShuffleSequences,
                                weight,
                                defaultTransparency,
                                defaultClosure
                            });
                        }
                    }
                    
                    */
                    // Figure 2 J-L
                    // Add combinations with different transparency values
                    for (double transparency : {0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0}) {
                        for (int repl = 0; repl < replications; ++repl) {
                            combinations.push_back({
                                adjMatrix, 
                                adjMatrixFlattened, 
                                strategy, 
                                defaultDistribution, 
                                defaultAlpha, 
                                repl, 
                                defaultSlope, 
                                defaultPayoffDist, 
                                usedShuffleSequences,
                                defaultEdgeWeight,
                                transparency,
                                defaultClosure
                            });
                        }
                    }
                    
                    
                }

            } else if (postfix == "maths" || postfix == "cook" || postfix == "honey" || postfix == "tuber") {
                double defaultSlope = (strategy == Strategy::Random || strategy == Strategy::Perfect) ? 0.0 : 2.0;

                auto transitiveClosure = computeTransitiveClosure(adjMatrix);
                auto transtiveClosureFlattened = adjMatrixToFlattenedString(transitiveClosure);
                
                for (int repl = 0; repl < replications; ++repl) {
                    combinations.push_back({
                        adjMatrix, 
                        adjMatrixFlattened, 
                        strategy, 
                        defaultDistribution, 
                        defaultAlpha, 
                        repl, 
                        defaultSlope, 
                        defaultPayoffDist, 
                        usedShuffleSequences,
                        defaultEdgeWeight,
                        defaultTransparency,
                        defaultClosure
                    });
                }
                
                for (int repl = 0; repl < replications; ++repl) {
                    combinations.push_back({
                        adjMatrix, 
                        adjMatrixFlattened, 
                        strategy, 
                        defaultDistribution, 
                        defaultAlpha, 
                        repl, 
                        defaultSlope, 
                        defaultPayoffDist, 
                        usedShuffleSequences,
                        defaultEdgeWeight,
                        10, //Transparency
                        defaultClosure
                    });
                }



                for (int repl = 0; repl < replications; ++repl) {
                    combinations.push_back({
                        transitiveClosure, 
                        transtiveClosureFlattened, 
                        strategy, 
                        defaultDistribution, 
                        defaultAlpha, 
                        repl, 
                        defaultSlope, 
                        defaultPayoffDist, 
                        usedShuffleSequences,
                        0.5, // Edge weight
                        defaultTransparency,
                        1 // Use transitive closure 
                    });
                }

            // Parameter settings used for figure 2 B-D 
            // Expected inputs: adj_mat_20.csv, adj_mat_30.csv, adj_mat_50.csv 
            } else {
                // For n > 8, use only the default slope
                double defaultSlope = (strategy == Strategy::Random || strategy == Strategy::Perfect) ? 0.0 : 2.0;
                
                for (int repl = 0; repl < replications; ++repl) {
                    combinations.push_back({
                        adjMatrix, 
                        adjMatrixFlattened, 
                        strategy, 
                        defaultDistribution, 
                        defaultAlpha, 
                        repl, 
                        defaultSlope, 
                        defaultPayoffDist, 
                        usedShuffleSequences,
                        defaultEdgeWeight,
                        defaultTransparency,
                        defaultClosure
                    });
                }
            }
        }
    }
    
    return combinations;
}

std::string stateToString(const Repertoire& state) {
    std::string result;
    result.reserve(state.size());
    for (double value : state) {
        result += value == 1.0 ? '1' : '0';
    }
    return result;
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

AdjacencyMatrix adjustMatrixWeights(const AdjacencyMatrix& adjMatrix, double edgeWeight) {
    auto weightedMatrix = adjMatrix;
    for (size_t row = 0; row < adjMatrix.size(); row++) {
        for (size_t col = 0; col < adjMatrix.size(); col++) {
            // if there's an edge, adjust the edge weight
            weightedMatrix[row][col] = adjMatrix[row][col] == 1.0 ? edgeWeight : 0.0;
        }
    }
    return weightedMatrix;
}

