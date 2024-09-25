#ifndef GRAPH_HPP 
#define GRAPH_HPP

#include "Types.hpp" 
#include <string>

std::vector<int> computeDistances(const AdjacencyMatrix& adjMatrix, Trait root);

std::vector<Trait> parentTraits(const AdjacencyMatrix& adjMatrix, Trait trait);

std::string adjMatrixToBinaryString(const AdjacencyMatrix& adjMatrix);

#endif // GRAPH_HPP