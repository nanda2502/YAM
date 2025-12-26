#include <atomic>
#include <iostream>
#include <numeric>
#include <vector>
#include <string>
#include <filesystem>

#include "Utils.hpp"
#include "Types.hpp"
#include "ModelConfig.hpp"
#include "IO.hpp"
#include "StringUtils.hpp"


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

    
        std::vector<ParamCombination> combinations = makeCombinations(adjacencyMatrices,replications, postfix);
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
    
        std::string csvHeader = "num_nodes,adj_mat,alpha,strategy,repl,steps,step_payoff,step_transitions,slope,distribution,absorbing,payoffdist,edge_weight,transparency,closure";
        std::vector<std::string> csvData;
        csvData.push_back(csvHeader);

        for (size_t i = 0; i < accumulatedResults.size(); ++i) {
            const AccumulatedResult& accumResult = accumulatedResults[i];
            if (accumResult.count > 0) {
                for (size_t step = 0; step < 20; step++) {
                const ParamCombination& comb = combinations[i];
                std::string formattedResult = formatResults(
                    comb.adjMatrix.size(),
                    adjMatrixToFlattenedString(comb.adjMatrix),
                    comb.alpha,
                    comb.strategy,
                    comb.repl,
                    step + 1,// add 1 since the index is 0-based
                    accumResult.totalExpectedPayoffPerStep[step],
                    accumResult.totalExpectedTransitionsPerStep[step],
                    comb.slope,
                    comb.distribution,
                    accumResult.absorbing,
                    comb.payoffDist,
                    comb.edgeWeight,
                    comb.transparency,
                    comb.closure
                );
                csvData.push_back(formattedResult);
                }

            }
        }

        writeCSV(outputDir, adj_idx, csvData);
       
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
    
    return 0;
}