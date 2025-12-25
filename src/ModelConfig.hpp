#ifndef MODELCONFIG_HPP
#define MODELCONFIG_HPP

#include <atomic>
#include <vector>
#include "Types.hpp"

// Model-specific configuration and parameter generation

// Generate parameter combinations for simulation runs
// This function defines the specific parameter sweeps used in the research article
std::vector<ParamCombination> makeCombinations(
    const std::vector<AdjacencyMatrix>& adjacencyMatrices, 
    int replications,
    const std::string& postfix
);

// Return slope values for a given strategy
std::vector<double> returnSlopeVector(Strategy strategy);

// Generate permutation sequences for payoff assignments
std::vector<std::vector<size_t>> makeShuffles(int n);

// Check if adjacency matrix is unconstrained (no prerequisites after root)
bool isUnconstrained(const AdjacencyMatrix& adjMatrix);

// Compute transitive closure of adjacency matrix using Floyd-Warshall
AdjacencyMatrix computeTransitiveClosure(const AdjacencyMatrix& adjacencyMatrix);

// Adjust edge weights in adjacency matrix
AdjacencyMatrix adjustMatrixWeights(const AdjacencyMatrix& adjMatrix, double edgeWeight);

// Process a single replication of the simulation
// Accumulates results across shuffle sequences for a given parameter combination
void processRepl(
    const ParamCombination& params,
    AccumulatedResult& accumResult,
    std::atomic<int>& failureCount
);

#endif // MODELCONFIG_HPP
