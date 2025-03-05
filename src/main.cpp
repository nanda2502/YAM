#include <algorithm>
#include <atomic>
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>

#include "Utils.hpp"
#include "Types.hpp"
#include "ExpectedSteps.hpp"


void processRepl(
    const AdjacencyMatrix& adjMatrix,
    const Strategy& strategy,
    double alpha,
    double slope,
    std::vector<std::vector<size_t>>& shuffleSequences,
    AccumulatedResult& accumResult,
    std::atomic<int>& failureCount,
    int adj_int,
    traitDistribution distribution
) {
    std::vector<double> totalExpectedPayoffPerStep(20, 0.0);
    std::vector<double> totalExpectedTransitionsPerStep(20, 0.0);
    std::vector<double> totalExpectedVariation(20, 0.0);
    double totalTimeToAbsorption = 0.0;
    int successCount = 0;

    // adj_int 0 is the unconstrained structure, which means that there is no need to shuffle the payoffs
    if (adj_int == 0) {
        shuffleSequences = {shuffleSequences[0]};
    }

    for (const auto& shuffleSequence : shuffleSequences) {
        double timeToAbsorption;
        timeToAbsorption = shuffleSequence == shuffleSequences[0] ? 0.0 : -1.0; // only compute time to absorption for the first shuffle sequence
        std::vector<double> expectedPayoffPerStep(20, 0.0);
        std::vector<double> expectedTransitionsPerStep(20, 0.0);
        std::vector<double> expectedVariation(20, 0.0);
        std::vector<std::vector<double>> transitionMatrix;

        if (computeExpectedSteps(adjMatrix, strategy, alpha, shuffleSequence,
                                 slope, expectedPayoffPerStep,
                                 expectedTransitionsPerStep, expectedVariation,
                                 transitionMatrix, distribution,
                                 timeToAbsorption)) {
          for (size_t i = 0; i < 20; ++i) {
            totalTimeToAbsorption += timeToAbsorption;
            totalExpectedPayoffPerStep[i] += expectedPayoffPerStep[i];
            totalExpectedTransitionsPerStep[i] += expectedTransitionsPerStep[i];
            totalExpectedVariation[i] += expectedVariation[i];
          }
          successCount++;
        } else {
          failureCount++;
        }
    }

    if (successCount > 0) {
        for (size_t i = 0; i < 20; ++i) {
            accumResult.count++;
            accumResult.totalExpectedPayoffPerStep[i] += totalExpectedPayoffPerStep[i] / shuffleSequences.size();
            accumResult.totalExpectedTransitionsPerStep[i] += totalExpectedTransitionsPerStep[i] / shuffleSequences.size();
            accumResult.totalExpectedVariation[i] += totalExpectedVariation[i] / shuffleSequences.size();
        }
        accumResult.absorbing += totalTimeToAbsorption / shuffleSequences.size();
    }
}

size_t factorial(size_t num) {
    size_t result = 1;
    for (size_t i = 2; i <= num; ++i) {
        result *= i;
    }
    return result;
}

int main(int argc, char* argv[]) {
    int n;
    int adj_int = parseArgs(argc, argv, n);
    int replications = 1;
    try {
        std::vector<double> alphas = {0.0};
        std::vector<Strategy> strategies = {
            Strategy::Random,
            Strategy::Payoff,
            Strategy::Proximal,
            Strategy::Prestige,
            Strategy::Conformity,
            Strategy::Perfect
        };

        std::vector<traitDistribution> distributions = {
            traitDistribution::Learnability//,
            //traitDistribution::Uniform,
            //traitDistribution::Depth,
            //traitDistribution::Shallowness,
            //traitDistribution::Payoffs
        };
    
        std::string outputDir = "../output";
        if (!std::filesystem::exists(outputDir)) {
            std::filesystem::create_directory(outputDir);
        }
    
        std::vector<AdjacencyMatrix> adjacencyMatricesAll = readAdjacencyMatrices(n);
        std::vector<AdjacencyMatrix> adjacencyMatrices(1, adjacencyMatricesAll[adj_int]);

        std::vector<size_t> perm(n - 1);
        std::iota(perm.begin(), perm.end(), 0);
        size_t sequenceCount = factorial(n - 1);
        std::vector<std::vector<size_t>> shuffleSequences;
        shuffleSequences.reserve(sequenceCount);  

        do {
            shuffleSequences.push_back(perm);
        } while (std::ranges::next_permutation(perm).found); 
    
    
        std::vector<ParamCombination> combinations = makeCombinations(adjacencyMatrices, strategies, alphas, replications, distributions);
        std::vector<AccumulatedResult> accumulatedResults(combinations.size());
        std::vector<std::atomic<int>> failureCounts(combinations.size());
    
        std::vector<size_t> indices(combinations.size());
        std::iota(indices.begin(), indices.end(), 0);
    
        #pragma omp parallel for
        for (unsigned long idx : indices) {
             const ParamCombination& comb = combinations[idx];
            processRepl(
                comb.adjMatrix,
                comb.strategy,
                comb.alpha,
                comb.slope,
                shuffleSequences,
                accumulatedResults[idx],
                failureCounts[idx],
                adj_int,
                comb.distribution
            );
        }
    

        std::string csvHeader = "num_nodes,adj_mat,alpha,strategy,repl,steps,step_payoff,step_transitions,step_variation,slope,distribution";
        std::vector<std::string> csvData;
        csvData.push_back(csvHeader);

        for (size_t i = 0; i < accumulatedResults.size(); ++i) {
            const AccumulatedResult& accumResult = accumulatedResults[i];
            if (accumResult.count > 0) {
                for (size_t step = 0; step < 20; step++) {
                const ParamCombination& comb = combinations[i];
                std::string formattedResult = formatResults(
                    n,
                    adjMatrixToBinaryString(comb.adjMatrix),
                    comb.alpha,
                    comb.strategy,
                    comb.repl,
                    step + 1,// add 1 since the index is 0-based
                    accumResult.totalExpectedPayoffPerStep[step],
                    accumResult.totalExpectedTransitionsPerStep[step],
                    accumResult.totalExpectedVariation[step],
                    comb.slope,
                    comb.distribution,
                    accumResult.absorbing
                );
                csvData.push_back(formattedResult);
                }

            }
        }

        writeAndCompressCSV(outputDir, adj_int, csvData);
       
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
    
    return 0;
}