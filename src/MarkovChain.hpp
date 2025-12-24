/*
#ifndef MARKOVCHAIN_HPP
#define MARKOVCHAIN_HPP

#include "Debug.hpp"
#include "Learning.hpp"
#include "Graph.hpp"
#include "Payoffs.hpp"
#include "LinAlg.hpp"
#include "Types.hpp"
#include "Utils.hpp"
#include "ExpectedSteps.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <numeric>
#include <cmath>

class MarkovChain {
public:
    // Constructor
    MarkovChain(
        const ParamCombination& params, 
        const std::vector<size_t>& shuffleSequence,                      
        std::vector<double>& expectedPayoffPerStep,
        std::vector<double>& expectedTransitionsPerStep,
        std::vector<double>& expectedVariation,                       
        std::vector<std::vector<double>>& transitionMatrix,
        double& timeToAbsorption,
        double& stationaryVariation
    );

    // Main computation method
    bool compute();

    // Public access to results (references to external vectors)
    std::vector<double>& expectedPayoffPerStep;
    std::vector<double>& expectedTransitionsPerStep;
    std::vector<double>& expectedVariation;
    double& timeToAbsorption;
    double& stationaryVariation;
    std::vector<std::vector<double>>& transitionMatrix;

private:
    // ========== ORIGINAL PARAMETERS ==========
    AdjacencyMatrix adjMatrix;
    Strategy strategy;
    double alpha;
    const std::vector<size_t>& shuffleSequence;
    double slope;
    double transparency;
    int payoffDist;
    traitDistribution distribution;

    // ========== DERIVED BASIC PROPERTIES ==========
    Strategy baseStrategy;
    size_t rootNode;
    size_t n;  // number of traits
    std::vector<double> distances;
    PayoffVector payoffs;
    Parents parents;
    std::vector<double> traitFrequencies;
    std::vector<Repertoire> allStates;

    // ========== REPERTOIRE MANAGEMENT ==========
    std::vector<Repertoire> repertoiresList;
    std::vector<Repertoire> finalRepertoiresList;
    std::unordered_map<Repertoire, int, RepertoireHash> repertoireIndexMap;
    std::unordered_map<Repertoire, int, RepertoireHash> finalRepertoireIndexMap;
    std::vector<std::pair<Repertoire, int>> repertoiresWithIndices;

    // ========== FREQUENCY DATA ==========
    std::unordered_map<Repertoire, double, RepertoireHash> stateFrequencies;
    std::unordered_map<Repertoire, double, RepertoireHash> initialStateFrequencies;
    std::unordered_map<Repertoire, double, RepertoireHash> inferredStateFrequencies;
    std::unordered_map<Repertoire, int, RepertoireHash> transientStateIndices;
    std::vector<Repertoire> transientStates;

    // ========== TRANSITION DATA ==========
    std::vector<Transition> allTransitions;
    std::vector<Transition> finalAllTransitions;

    // ========== MATRIX COMPUTATION STATE ==========
    std::vector<std::vector<double>> preliminaryTransitionMatrix;
    std::vector<std::vector<double>> reorderedTransitionMatrix;
    std::vector<std::vector<double>> fundamentalMatrix;
    std::vector<std::vector<double>> iMinusQ;
    std::vector<std::vector<double>> LU;
    std::vector<int> p;  // permutation vector for LU decomposition
    std::unordered_map<int, int> oldToNewIndexMap;
    int numTransientStates;

    // ========== PAYOFF INFORMATION ==========
    std::vector<double> statePayoffs;
    std::vector<double> allStatesPayoffs;
    std::vector<double> initialStatePayoffs;

    // ========== KEY REFERENCES ==========
    int initialStateIndex;
    Repertoire initialRepertoire;
    Repertoire absorbingState;

    // ========== PRIVATE MEMBER FUNCTIONS ==========
    
    // Main computation phases
    void buildInitialTransitionMatrix();
    void computeFundamentalMatrix();
    void calculateStateFrequencies();
    void updateTraitFrequencies();
    void buildFinalTransitionMatrix();
    void computeResults();

    // Debug utilities
    void debugPrint() const;
};
#endif // MARKOVCHAIN_HPP
*/