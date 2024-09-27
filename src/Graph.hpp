#ifndef GRAPH_HPP 
#define GRAPH_HPP

#include "Types.hpp" 

std::vector<int> computeDistances(const AdjacencyMatrix& adjMatrix, Trait root);

std::vector<Trait> parentTraits(const AdjacencyMatrix& adjMatrix, Trait trait);

#endif // GRAPH_HPP