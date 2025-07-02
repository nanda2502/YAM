#include "Graph.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <queue>

std::vector<double> computeDistances(const AdjacencyMatrix& adjMatrix, Trait root) {
    size_t n = adjMatrix.size();
    std::vector<double> distances(n, 0.0);
    
    // Check if any edge weight is not 0.0 or 1.0 (weighted graph)
    bool isWeighted = std::any_of(adjMatrix.begin(), adjMatrix.end(),
        [](const auto& row) {
            return std::any_of(row.begin(), row.end(),
                [](double w) { return w > 0.0 && w != 1.0; });
        });
    
    if (isWeighted) {
        // Weighted: sum of incoming edge weights
        for (size_t trait = 0; trait < n; ++trait) {
            for (size_t source = 0; source < n; ++source) {
                distances[trait] += adjMatrix[source][trait];
            }
        }
    } else {
        // Unweighted: BFS from root (your existing logic)
        std::vector<int> intDistances(n, -1);
        std::queue<Trait> q;
        
        intDistances[root] = 0;
        q.push(root);
        
        while (!q.empty()) {
            Trait current = q.front();
            q.pop();
            
            for (size_t neighbor = 0; neighbor < n; ++neighbor) {
                if (adjMatrix[current][neighbor] > 0.0 && intDistances[neighbor] == -1) {
                    intDistances[neighbor] = intDistances[current] + 1;
                    q.push(neighbor);
                }
            }
        }
        
        // Convert to double and check connectivity
        for (size_t node = 0; node < n; ++node) {
            if (intDistances[node] == -1) {
                throw std::runtime_error("Graph is not connected");
            }
            distances[node] = static_cast<double>(intDistances[node]);
        }
    }
    
    return distances;
}

std::vector<Trait> parentTraits(const AdjacencyMatrix& adjMatrix, Trait trait) {
    std::vector<Trait> parents;
    size_t n = adjMatrix.size();

    for (size_t i = 0; i < n; ++i) {
        if (adjMatrix[i][trait] == 1.0) {
            parents.push_back(i);
        }
    }
    return parents;
}

