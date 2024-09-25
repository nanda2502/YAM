#ifndef EXPECTEDSTEPS_HPP
#define EXPECTEDSTEPS_HPP

#include "Types.hpp"
#include "Learning.hpp"

#include <random>
#include <vector>
#include <unordered_map>
#include <tuple>


std::pair<double, std::vector<std::vector<double>>> computeExpectedSteps(
    const AdjacencyMatrix& adjacencyMatrix,
    Strategy strategy,
    double alpha,
    std::mt19937& gen
);

std::vector<std::vector<double>> buildTransitionMatrix(
    const std::vector<Repertoire>& repertoiresList,
    const std::unordered_map<Repertoire, int, RepertoireHash>& repertoireIndexMap,
    Strategy strategy,
    const AdjacencyMatrix& adjacencyMatrix,
    const PayoffVector& payoffs
);

std::tuple<std::vector<std::vector<double>>, std::unordered_map<int, int>, int> reorderTransitionMatrix(
    const std::vector<std::vector<double>>& transitionMatrix,
    const std::vector<std::pair<Repertoire, int>>& repertoiresWithIndices,
    const std::unordered_map<Repertoire, int, RepertoireHash>& repertoireIndexMap,
    Trait rootNode
);

double computeExpectedStepsFromMatrix(
    const std::vector<std::vector<double>>& reorderedTransitionMatrix,
    int numTransientStates,
    int initialStateNewIndex
);

#endif // EXPECTEDSTEPS_HPP