#include "Graph.hpp"

#include <cstddef>
#include <stdexcept>
#include <queue>

std::vector<int> computeDistances(const AdjacencyMatrix& adjMatrix, Trait root) {
    size_t n = adjMatrix.size();
    std::vector<int> distances(n, -1);
    std::queue<Trait> q;

    distances[root] = 0;
    q.push(root);

    while (!q.empty()) {
        Trait current = q.front();
        q.pop();

        for (size_t neighbor = 0; neighbor < n; ++neighbor) {
            if (adjMatrix[current][neighbor] == 1.0 && distances[neighbor] == -1) {
                distances[neighbor] = distances[current] + 1;
                q.push(neighbor);
            }
        }
    }

    for (size_t node = 0; node < n; ++node) {
        if (distances[node] == -1) {
            throw std::runtime_error("Graph is not connected");
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

