#include "Learning.hpp"
#include "Debug.hpp"
#include "Types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <queue>
#include <unordered_map>
#include <unordered_set>

std::vector<double> normalizedMissingPrereqs(
    const Repertoire& repertoire,
    const AdjacencyMatrix& adjMatrix
) {
    size_t n = repertoire.size();
    std::vector<double> ancestorCounts(n, 0.0);
    
    // For each unlearned trait, count unknown ancestors
    for (size_t trait = 1; trait < n; ++trait) {
        if (repertoire[trait] == 1.0) {
            continue; // Skip learned traits
        }
        
        // BFS/DFS to find all ancestors of this trait
        std::queue<size_t> queue;
        std::unordered_set<size_t> visited;
        std::unordered_set<size_t> ancestors;
        
        queue.push(trait);
        visited.insert(trait);
        
        while (!queue.empty()) {
            size_t current = queue.front();
            queue.pop();
            
            // Find all parents of current node
            for (size_t parent = 0; parent < n; ++parent) {
                if (adjMatrix[parent][current] > 0.0) {
                    ancestors.insert(parent);
                    
                    if (visited.find(parent) == visited.end()) {
                        visited.insert(parent);
                        queue.push(parent);
                    }
                }
            }
        }
        
        // Count unknown ancestors
        for (size_t ancestor : ancestors) {
            if (repertoire[ancestor] == 0.0) {
                ancestorCounts[trait] += 1.0; // or could weight by adjMatrix[ancestor][...] if desired
            }
        }
        
        DEBUG_PRINT(2, "Trait " << trait << " has " << ancestorCounts[trait] << " unknown ancestors");
    }
    
    // Normalize
    double total = std::accumulate(ancestorCounts.begin(), ancestorCounts.end(), 0.0);
    if (total > 0.0) {
        std::ranges::transform(ancestorCounts, ancestorCounts.begin(),
            [total](double w) { return w / total; });
    }
    
    return ancestorCounts;
}


std::vector<double> learnability(
    const Repertoire& repertoire,
    const AdjacencyMatrix& adjMatrix
) {
    size_t n = repertoire.size();
    std::vector<double> learnable(n, 0.0);
    
    for (Trait trait = 0; trait < n; ++trait) {
        // If trait is already learned, it's not learnable
        if (repertoire[trait] == 1.0) {
            learnable[trait] = 0.0;
            continue;
        }
        
        // Start with probability 1.0
        double probability = 1.0;
        
        // Check all potential parents
        for (Trait parent = 0; parent < n; ++parent) {
            double edgeWeight = adjMatrix[parent][trait];
            
            if (edgeWeight > 0.0) {  // This is a direct parent
                if (repertoire[parent] == 0.0) {
                    // Parent is unknown - contribute with probability (1-edgeWeight)
                    // This represents the chance that this unknown parent isn't required
                    probability *= (1.0 - edgeWeight);
                }
                // If parent is known, it doesn't reduce the probability (contributes with factor 1.0)
            }
        }
                
        learnable[trait] = probability;
    }
    
    return learnable;
}

double computeDelta(const Repertoire& r, const Repertoire& s) {
    // count the number of traits that are present in target state s but not in current state r
    return std::inner_product(s.begin(), s.end(), r.begin(), 0.0,
        std::plus<>(), [](double s_i, double r_i) { return (s_i == 1.0 && r_i == 0.0) ? 1 : 0; });
}

std::vector<double> proximalBaseWeights(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope
) {
    // Step 1: Calculate raw bias for each state
    std::vector<double> stateBiases(allStates.size(), 0.0);
    std::vector<bool> validStates(allStates.size(), false);
    
    for (size_t i = 0; i < allStates.size(); ++i) {
        const auto& state = allStates[i];
        double delta = computeDelta(repertoire, state);
        if (delta > 0) {
            stateBiases[i] = std::pow(delta, -slope);
            validStates[i] = true;
        }
    }
    
    // Step 2: Normalize state biases
    double totalBias = 0.0;
    int validCount = 0;
    for (size_t i = 0; i < allStates.size(); ++i) {
        if (validStates[i]) {
            totalBias += stateBiases[i];
            validCount++;
        }
    }

    double meanBias = totalBias / validCount;
    if (meanBias > 0.0) {
        for (double& bias : stateBiases) {
            bias /= meanBias;  // Now average bias = 1.0
        }
    }
    
    // Step 3: Apply state frequencies to get trait weights
    std::vector<double> traitWeights(repertoire.size(), 0.0);
    for (size_t i = 0; i < allStates.size(); ++i) {
        if (!validStates[i]) continue;
        
        // Get the frequency of the current state
        const auto& state = allStates[i];
        auto it = stateFrequencies.find(state);
        if (it == stateFrequencies.end()) continue;
        double stateFreq = it->second;

        double normalizedBias = stateBiases[i];
        
        for (Trait trait = 0; trait < repertoire.size(); ++trait) {
            if (repertoire[trait] == 0.0 && state[trait] == 1.0) {
                traitWeights[trait] += normalizedBias * stateFreq;
            }
        }
    }
    
    return traitWeights;
}

std::vector<double> payoffBaseWeights(
    const std::vector<double>& payoffs, 
    const std::vector<double>& traitFrequencies, 
    double slope
) {
    // Step 1: Calculate raw biases  
    std::vector<double> rawBiases(payoffs.size());
    for (size_t i = 0; i < payoffs.size(); ++i) {
        rawBiases[i] = std::pow(payoffs[i], slope);
    }
    
    // Step 2: Normalize biases to mean = 1.0 (exclude zero-payoff traits)
    double totalBias = 0.0;
    int validCount = 0;
    for (size_t i = 0; i < payoffs.size(); ++i) {
        if (payoffs[i] > 0.0) {  // Only count non-zero payoffs
            totalBias += rawBiases[i];
            validCount++;
        }
    }
    
    double meanBias = (validCount > 0) ? totalBias / validCount : 1.0;
    
    // Step 3: Apply normalization and frequencies
    std::vector<double> finalWeights(payoffs.size());
    for (size_t i = 0; i < payoffs.size(); ++i) {
        double normalizedBias = rawBiases[i] / meanBias;
        finalWeights[i] = normalizedBias * traitFrequencies[i];
    }
    
    return finalWeights;
}

std::vector<double> prestigeBaseWeights(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    const std::vector<double>& statePayoffs,
    double slope
) {
    // Step 1: Calculate raw bias for each state
    std::vector<double> stateBiases(allStates.size(), 0.0);
    std::vector<bool> validStates(allStates.size(), false);
    
    for (size_t i = 0; i < allStates.size(); ++i) {
        const auto& state = allStates[i];
        double delta = computeDelta(repertoire, state);
        if (delta > 0) {
            stateBiases[i] = std::pow(statePayoffs[i], slope);
            validStates[i] = true;
        }
    }
    
    // Step 2: Normalize state biases
    double totalBias = 0.0;
    int validCount = 0;
    for (size_t i = 0; i < allStates.size(); ++i) {
        if (validStates[i]) {
            totalBias += stateBiases[i];
            validCount++;
        }
    }

    double meanBias = totalBias / validCount;
    if (meanBias > 0.0) {
        for (double& bias : stateBiases) {
            bias /= meanBias;  // Now average bias = 1.0
        }
    }
    
    // Step 3: Apply state frequencies to get trait weights
    std::vector<double> traitWeights(repertoire.size(), 0.0);
    for (size_t i = 0; i < allStates.size(); ++i) {
        if (!validStates[i]) continue;
        
        const auto& state = allStates[i];
        auto it = stateFrequencies.find(state);
        if (it == stateFrequencies.end()) continue;
        
        double normalizedBias = stateBiases[i];
        double stateFreq = it->second;
        
        for (Trait trait = 0; trait < repertoire.size(); ++trait) {
            if (repertoire[trait] == 0.0 && state[trait] == 1.0) {
                traitWeights[trait] += normalizedBias * stateFreq;
            }
        }
    }
    
    return traitWeights;
}

std::vector<double> conformityBaseWeights(
    const std::vector<double>& traitFrequencies,
    double slope
) {
    std::vector<double> w_star(traitFrequencies.size());
    
    // Conformity: frequency^slope 
    std::ranges::transform(traitFrequencies, w_star.begin(), 
        [slope](double f) {
            return std::pow(f, slope);
        });
    
    return w_star;
}

std::vector<double> anticonformityBaseWeights(
    const std::vector<double>& traitFrequencies,
    double slope
) {
    
    std::vector<double> w_star(traitFrequencies.size());
    std::ranges::transform(traitFrequencies, w_star.begin(),
        [slope](double f) {
            return 1.0 - std::pow(f, slope); // 1 - frequency^slope
        });
    return w_star;
}

std::vector<double> perfectBaseWeights(
    const Repertoire& repertoire,
    const PayoffVector& payoffs,  
    const std::vector<double>& learnableProbs
) {
    
    std::vector<double> w_star(repertoire.size(), 0.0);
    // Always learn the trait with the highest expected payoff that is not already learned
    double maxExpectedPayoff = 0.0;
    size_t bestTrait = 0;
    for (size_t trait = 0; trait < repertoire.size(); ++trait) {
        double expectedPayoff = payoffs[trait] * learnableProbs[trait];
        if (expectedPayoff > maxExpectedPayoff) {
            maxExpectedPayoff = expectedPayoff;
            bestTrait = trait;
        }
    }
    // Set the weight for the best trait to 1.0, others to 0.0
    w_star[bestTrait] = 1.0;
    

    return w_star;
}

std::vector<double> baseWeights(
    Strategy strategy,
    const Repertoire& repertoire,
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope,
    const std::vector<double>& statePayoffs,
    const std::vector<double>& learnableProbs   
) {
    switch (strategy) {
    case Random:
        return traitFrequencies;
    case Payoff:
        return payoffBaseWeights(payoffs, traitFrequencies, slope);
    case Proximal:
        return proximalBaseWeights(repertoire, stateFrequencies, allStates, slope);
    case Prestige:
        return prestigeBaseWeights(repertoire, stateFrequencies, allStates, statePayoffs, slope);
    case Conformity:
        return conformityBaseWeights(traitFrequencies, slope);
    case Anticonformity:
        return anticonformityBaseWeights(traitFrequencies, slope);
    case Perfect:
        return perfectBaseWeights(repertoire, payoffs, learnableProbs);
    default:
        throw std::runtime_error("Unknown strategy");
    }
}

std::vector<double> normalizedWeights(
    Strategy strategy,
    const Repertoire& repertoire,
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope,
    double lambda,
    const std::vector<double>& statePayoffs,
    const std::vector<double>& learnableProbs,
    const AdjacencyMatrix& adjMatrix
)  {
    std::vector<double> w_star = baseWeights(strategy, repertoire, payoffs, traitFrequencies, stateFrequencies, allStates, slope, statePayoffs, learnableProbs);

    DEBUG_PRINT(2, "Current repertoire:");
    if (DEBUG_LEVEL >= 2) {
        for (double i: repertoire) {
            std::cout << (i == 1.0 ? "1" : "0");
        }
        std::cout << '\n';
    }

    DEBUG_PRINT(2, "Base weights:");
    if (DEBUG_LEVEL >= 2) {
        for (size_t i = 0; i < w_star.size(); ++i) {
            std::cout << "Trait " << i <<": " << w_star[i] << '\n';
        }
    }

    std::vector<double> w_unlearned(repertoire.size());
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        w_unlearned[trait] = (repertoire[trait] == 1.0 || traitFrequencies[trait] == 0.0) ? 0.0 : w_star[trait];
    }

    DEBUG_PRINT (2, "Weights after setting learned ones to 0:");
    if (DEBUG_LEVEL >= 2) {
        for (size_t i = 0; i < w_unlearned.size(); ++i) {
            std::cout << "Trait " << i <<": " << w_unlearned[i] << '\n';
        }
    }

    // Attempt weights depend on partial structure knowledge. When lambda is 0, structure is completely opaque. 
    auto missingPrereqs = normalizedMissingPrereqs(repertoire, adjMatrix);

    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        w_unlearned[trait] = w_unlearned[trait] * std::pow((1 - missingPrereqs[trait]), lambda);
    }
    
    double total = std::accumulate(w_unlearned.begin(), w_unlearned.end(), 0.0);

    if (total == 0.0) {
        return std::vector<double>(repertoire.size(), 0.0);
    }

    std::ranges::transform(w_unlearned, w_unlearned.begin(),
        [total](double w) { return w / total; });

    DEBUG_PRINT(2, "Normalized weights:");
    if (DEBUG_LEVEL >= 2) {
        for (size_t i = 0; i < w_unlearned.size(); ++i) {
            std::cout << "Trait " << i <<": " << w_unlearned[i] << '\n';
        }
    };

    if (strategy == Strategy::Prestige && DEBUG_LEVEL >= 1) {
        std::cout << "After normalization - trait probabilities:\n";
        for (size_t i = 1; i < w_unlearned.size(); ++i) {
            std::cout << "  Trait " << i << ": " << w_unlearned[i] << std::endl;
        }
    }

    return w_unlearned;
}

Repertoire learnTrait(const Repertoire& repertoire, Trait trait) {
    Repertoire newRepertoire = repertoire;
    newRepertoire[trait] = 1.0;
    return newRepertoire;
}

std::vector<std::pair<Repertoire, double>> transitionFromState(
    Strategy strategy,
    const Repertoire& repertoire, 
    const PayoffVector& payoffs, 
    const std::vector<double>& traitFrequencies,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope,
    double lambda,
    const AdjacencyMatrix& adjMatrix,
    const std::vector<double>& statePayoffs
) {
    std::vector<Repertoire> newStates = retrieveBetterRepertoires(allStates, repertoire);
    std::vector<double> learnableProbs = learnability(repertoire, adjMatrix);
    std::vector<double> w = normalizedWeights(strategy, repertoire, payoffs, traitFrequencies, stateFrequencies, newStates, slope, lambda, statePayoffs, learnableProbs, adjMatrix);
    

    std::vector<std::pair<Repertoire, double>> transitions;

    for(Trait trait = 0; trait < repertoire.size(); ++trait) {
        // Multiply the normalized weight by the probability that the trait is learnable
        double transitionProb = w[trait] * learnableProbs[trait];
        if (transitionProb > 0.0) {
            transitions.emplace_back(learnTrait(repertoire, trait), transitionProb);
        }
    }

    return transitions;
}

double stayProbability(std::vector<std::pair<Repertoire, double>> transitions) {

    double totalTransitionProbability = std::accumulate(
        transitions.begin(), transitions.end(), 0.0,
        [](double sum, const auto& transition) { 
            return sum + transition.second;
        }
    );

    return 1.0 - totalTransitionProbability;
}

std::pair<std::vector<Repertoire>, std::vector<std::vector<std::pair<Repertoire, double>>>>  generateReachableRepertoires(
    Strategy strategy, 
    const AdjacencyMatrix& adjMatrix, 
    const PayoffVector& payoffs, 
    const std::vector<double>& traitFrequencies,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope,
    const std::vector<double>& statePayoffs,
    double lambda
) {
    size_t n = adjMatrix.size();
    Repertoire initialRepertoire(n, 0.0);
    initialRepertoire[0] = 1.0; // root trait is always learned

    std::queue<Repertoire> queue;
    std::unordered_set<Repertoire, RepertoireHash> visited;
    std::vector<Repertoire> result;
    std::vector<std::vector<std::pair<Repertoire, double>>> allTransitions;

    queue.push(initialRepertoire);

    while (!queue.empty()) {
        Repertoire r = queue.front();
        queue.pop();

        if (visited.find(r) == visited.end()) {
            visited.insert(r);
            result.push_back(r);

            auto transitions = transitionFromState(strategy, r, payoffs, traitFrequencies, stateFrequencies, allStates, slope, lambda, adjMatrix, statePayoffs);
            allTransitions.push_back(transitions);

            for (const auto& transition : transitions) {
                const Repertoire& r_prime = transition.first;
                if (visited.find(r_prime) == visited.end()) {
                    queue.push(r_prime);
                }
            }
        }
    }

    return {result, allTransitions};
}

std::vector<Repertoire> generateAllRepertoires(const AdjacencyMatrix& adjMatrix) {
    size_t n = adjMatrix.size();
    Repertoire initialRepertoire(n, 0.0);
    initialRepertoire[0] = 1.0; // root trait is always learned

    std::queue<Repertoire> queue;
    std::unordered_set<Repertoire, RepertoireHash> visited;
    std::vector<Repertoire> result;

    queue.push(initialRepertoire);

    while (!queue.empty()) {
        Repertoire r = queue.front();
        queue.pop();

        if (visited.find(r) == visited.end()) {
            visited.insert(r);
            result.push_back(r);

            std::vector<double> learnableProbs = learnability(r, adjMatrix);
            for (Trait trait = 0; trait < n; ++trait) {
                if (learnableProbs[trait] > 0.0) {
                    Repertoire r_new = learnTrait(r, trait);
                    if (visited.find(r_new) == visited.end()) {
                        queue.push(r_new);
                    }
                }
            }
        }
    }
    return result;
}

std::vector<Repertoire> retrieveBetterRepertoires(const std::vector<Repertoire>& repertoires, const Repertoire& singleRepertoire) {
    std::vector<Repertoire> result;
    for (const Repertoire& r : repertoires) {
        for (size_t trait = 0; trait < r.size(); ++trait) {
            if (r[trait] == 1.0 && singleRepertoire[trait] == 0.0) {
                result.push_back(r);
                break;
            }
        }
    }
    return result;
}