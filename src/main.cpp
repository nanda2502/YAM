#include <atomic>
#include <iostream>
#include <numeric>
#include <vector>
#include <string>
#include <filesystem>

#include "Utils.hpp"
#include "Types.hpp"
#include "ExpectedSteps.hpp"


void processRepl(
    const ParamCombination& params,
    AccumulatedResult& accumResult,
    std::atomic<int>& failureCount
) {
    std::cout << "Processing strategy: " << strategyToString(params.strategy) << std::endl;
    std::vector<double> totalExpectedPayoffPerStep(20, 0.0);
    std::vector<double> totalExpectedTransitionsPerStep(20, 0.0);
    std::vector<double> totalExpectedVariation(20, 0.0);
    double totalTimeToAbsorption = 0.0;
    double totalStationaryVariation = 0.0;
    int successCount = 0;

    AdjacencyMatrix weightedMatrix;
    if (params.edgeWeight != 1.0) {
        weightedMatrix = adjustMatrix(params.adjMatrix, params.edgeWeight);
    } else {
        weightedMatrix = params.adjMatrix;
    }

    for (const auto& shuffleSequence : params.shuffleSequences) {
        double timeToAbsorption;
        double stationaryVariation;
        timeToAbsorption = shuffleSequence == params.shuffleSequences[0] ? 0.0 : -1.0; // only compute time to absorption for the first shuffle sequence
        std::vector<double> expectedPayoffPerStep(20, 0.0);
        std::vector<double> expectedTransitionsPerStep(20, 0.0);
        std::vector<double> expectedVariation(20, 0.0);
        std::vector<std::vector<double>> transitionMatrix;

        if (computeExpectedSteps(weightedMatrix, params.strategy, params.alpha, shuffleSequence,
                                 params.slope, params.lambda, params.payoffDist, params.distribution, 
                                 expectedPayoffPerStep,
                                 expectedTransitionsPerStep,
                                 expectedVariation,
                                 transitionMatrix, 
                                 timeToAbsorption,
                                 stationaryVariation)) {
            totalTimeToAbsorption += timeToAbsorption;
            totalStationaryVariation += stationaryVariation;
            for (size_t i = 0; i < 20; ++i) {
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
            accumResult.totalExpectedPayoffPerStep[i] += totalExpectedPayoffPerStep[i] / params.shuffleSequences.size();
            accumResult.totalExpectedTransitionsPerStep[i] += totalExpectedTransitionsPerStep[i] / params.shuffleSequences.size();
            accumResult.totalExpectedVariation[i] += totalExpectedVariation[i] / params.shuffleSequences.size();
        }
        accumResult.absorbing += totalTimeToAbsorption / params.shuffleSequences.size();
        accumResult.stationaryVariation += totalStationaryVariation / params.shuffleSequences.size();
    }

    std::cout << "Failures: " << failureCount.load() << std::endl;
}

int main(int argc, char* argv[]) {
    std::string postfix = "8"; // Default value
    int adj_idx = parseArgs(argc, argv, postfix);
    int replications = 1;
    try {
   
        std::string outputDir = "../output";
        if (!std::filesystem::exists(outputDir)) {
            std::filesystem::create_directory(outputDir);
        }
    
        std::vector<AdjacencyMatrix> adjacencyMatricesAll = readAdjacencyMatrices(postfix);
        std::vector<AdjacencyMatrix> adjacencyMatrices(1, adjacencyMatricesAll[adj_idx]);

    
        std::vector<ParamCombination> combinations = makeCombinations(adjacencyMatrices,replications);
        std::vector<AccumulatedResult> accumulatedResults(combinations.size());
        std::vector<std::atomic<int>> failureCounts(combinations.size());
    
        std::vector<size_t> indices(combinations.size());
        std::iota(indices.begin(), indices.end(), 0);
    
        #pragma omp parallel for
        for (unsigned long idx : indices) {
            processRepl(
                combinations[idx],
                accumulatedResults[idx],
                failureCounts[idx]
            );
        }
    
        std::string csvHeader = "num_nodes,adj_mat,alpha,strategy,repl,steps,step_payoff,step_transitions,step_variation,slope,distribution,absorbing,stationary_variation,payoffdist,edge_weight,lambda";
        std::vector<std::string> csvData;
        csvData.push_back(csvHeader);

        for (size_t i = 0; i < accumulatedResults.size(); ++i) {
            const AccumulatedResult& accumResult = accumulatedResults[i];
            if (accumResult.count > 0) {
                for (size_t step = 0; step < 20; step++) {
                const ParamCombination& comb = combinations[i];
                std::string formattedResult = formatResults(
                    comb.adjMatrix.size(),
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
                    accumResult.absorbing,
                    accumResult.stationaryVariation,
                    comb.payoffDist,
                    comb.edgeWeight,
                    comb.lambda
                );
                csvData.push_back(formattedResult);
                }

            }
        }

        writeAndCompressCSV(outputDir, adj_idx, csvData);
       
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
    
    return 0;
}