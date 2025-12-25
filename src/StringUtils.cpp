#include "StringUtils.hpp"
#include <iomanip>
#include <sstream>
#include <stdexcept>

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

std::string formatAdjMat(const std::string& adj_string, int n) {
    std::string adj_mat;
    int col_idx = 0;
    for (char i : adj_string) {
        if (col_idx == n) {
            adj_mat += '\n';
            col_idx = 0;
        }
        adj_mat += i;
        col_idx++;
    }
    return adj_mat;
}

std::string stateToString(const Repertoire& state) {
    std::string result;
    result.reserve(state.size());
    for (double value : state) {
        result += value == 1.0 ? '1' : '0';
    }
    return result;
}
