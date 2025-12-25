// This file implements the learning-related functions. 
// This includes the definitions of social learning strategies and the mechanism for state transitions between repertoires.

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
    const AdjacencyMatrix& adjMatrix,
    double transparency 
) {
    // Check if adjMatrix is strictly binary (0.0 or 1.0)
    for (size_t i = 0; i < adjMatrix.size(); ++i) {
        for (size_t j = 0; j < adjMatrix[i].size(); ++j) {
        double val = adjMatrix[i][j];
        if (val != 0.0 && val != 1.0 && transparency > 0.0) {
            throw std::runtime_error("Weighted matrices are not supported with transparency: adjMatrix[" +
                        std::to_string(i) + "][" + std::to_string(j) + "] = " + std::to_string(val));
        }
        }
    }

    size_t n = repertoire.size();
    std::vector<double> ancestorCounts(n, 0.0);
    
    // For each unlearned trait, count unknown ancestors
    for (size_t trait = 1; trait < n; ++trait) {
        if (repertoire[trait] == 1.0) {
            continue; // Skip learned traits
        }
        
        // BFS to find all ancestors of this trait
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
                ancestorCounts[trait] += 1.0; 
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
    const AdjacencyMatrix& adjMatrix // w[prereq][trait]
) {
    size_t n = repertoire.size();
    std::vector<double> learnable(n, 0.0);
    
    for (Trait trait = 0; trait < n; ++trait) {
         // Already known: not a candidate for learning attempts.
        if (repertoire[trait] == 1.0) {
            learnable[trait] = 0.0;
            continue;
        }
        
        // Learnability under soft dependencies: Π over missing prerequisites (1 - w).
        double probability = 1.0;
        
        // Check all potential prerequisites
        for (Trait prereq = 0; prereq < n; ++prereq) {
            double edgeWeight = adjMatrix[prereq][trait];
            
            if (edgeWeight > 0.0) {  // This is a direct prerequisite
                if (repertoire[prereq] == 0.0) {
                    // Prerequisite is not in repertoire - contribute with probability (1-edgeWeight) to learnability
                    // This represents the chance that this unknown parent isn't required
                    probability *= (1.0 - edgeWeight);
                }
                // If prerequisite is in repertoire, it doesn't reduce the probability (contributes with factor 1.0)
            }
        }
                
        learnable[trait] = probability;
    }
    
    return learnable;
}

// Count traits present in demo but missing in learner
int computeDelta(const Repertoire& learner, const Repertoire& demo) {
    int gap = 0;
    for (size_t t = 0; t < learner.size(); ++t) {
        if (demo[t] == 1.0 && learner[t] == 0.0) ++gap;
    }
    return gap;
}

std::vector<double> proximalBaseWeights(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    const std::vector<double>& payoffs,
    const std::vector<double>& depths,
    TraitDistribution dist,
    double slope
) {
    const size_t n = repertoire.size();
    std::vector<double> weights(n, 0.0);
    
    // Default case: traits are expressed with uniform probability
    if (dist == TraitDistribution::Uniform) {
        // Loop over each trait
        for (size_t trait = 0; trait < n; ++trait) {
            if (repertoire[trait] != 0) continue; // learner already has this trait

            // Loop over all potential demonstrator repertoires in the stateFrequencies map
            for (const auto& [demoRep, f_demo] : stateFrequencies) {
                if (demoRep[trait] == 0) continue; // demonstrator doesn't have the trait

                const auto delta = computeDelta(repertoire, demoRep);
                if (delta == 0) continue; // don't learn from demonstrators with no new traits
                
                // Proximal weight contribution: f_demo * delta^{-alpha}
                weights[trait] += f_demo * std::pow(delta, -slope);
        
            }
            
        }
        return weights;
    }

    auto sq = [](double x) { return x * x; };
    
    for (const auto& [demoRep, f_demo] : stateFrequencies) {
        const auto delta = computeDelta(repertoire, demoRep);
        if (delta == 0) continue;

        const double base = f_demo * std::pow(delta, -slope);

        // Normalize within this demonstrator so the average expression multiplier is ~1
        double sumSq = 0.0;
        size_t k = 0;
        for (size_t t = 0; t < n; ++t) {
            if (demoRep[t] == 0) continue;
            const double v = (dist == TraitDistribution::Depth) ? depths[t] : payoffs[t];
            sumSq += sq(v);
            ++k;
        }
        const double norm = (k > 0 && sumSq > 0.0) ? (sumSq / static_cast<double>(k)) : 1.0;

        for (size_t trait = 0; trait < n; ++trait) {
            if (repertoire[trait] != 0) continue;
            if (demoRep[trait] == 0) continue;

            // Expression bias: scale by trait attribute and normalize within demonstrator
            const double v = (dist == TraitDistribution::Depth) ? depths[trait] : payoffs[trait];
            const double expr = sq(v) / norm;

            weights[trait] += base * expr;
        }
    }

    return weights;
}



std::vector<double> payoffBaseWeights(
    const std::vector<double>& payoffs, 
    const std::vector<double>& traitFrequencies, 
    double slope
) {
    std::vector<double> weights(payoffs.size());
    
    for (size_t i = 0; i < payoffs.size(); ++i) {
        weights[i] = traitFrequencies[i] * std::pow(payoffs[i], slope);
    }
    
    return weights;
}

std::vector<double> prestigeBaseWeights(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    const std::vector<double>& statePayoffs,
    const std::vector<double>& payoffs,
    const std::vector<double>& depths,
    TraitDistribution dist,
    double slope
) {
    std::vector<double> weights(repertoire.size(), 0.0);

    // Map demonstrator repertoire -> total payoff (prestige signal).
    std::unordered_map<Repertoire, double, RepertoireHash> payoffByState;
    if (allStates.size() == statePayoffs.size()) {
        payoffByState.reserve(allStates.size());
        for (size_t i = 0; i < allStates.size(); ++i) {
            payoffByState.emplace(allStates[i], statePayoffs[i]);
        }
    }

    const size_t n = repertoire.size();

    // Default case: traits are expressed with uniform probability
    if (dist == TraitDistribution::Uniform) {
        for (size_t trait = 0; trait < n; ++trait) {
            if (repertoire[trait] != 0.0) continue; // learner already has this trait

            for (const auto& [demoState, f_demo] : stateFrequencies) {
                if (demoState[trait] != 1.0) continue; // demonstrator doesn't express the trait

                double demoPayoff = 0.0;
                if (!payoffByState.empty()) {
                    auto it = payoffByState.find(demoState);
                    if (it != payoffByState.end()) demoPayoff = it->second;
                } else {
                    for (size_t j = 0; j < demoState.size(); ++j) {
                        if (demoState[j] == 1.0) demoPayoff += payoffs[j];
                    }
                }

                const double prestigeBias = (slope == 0.0 ? 1.0 : std::pow(demoPayoff, slope));
                // Prestige weight: sum_{demo contains trait} f_demo * (totalPayoff_demo ^ slope)
                weights[trait] += f_demo * prestigeBias;
            }
        }

        return weights;
    }

    auto sq = [](double x) { return x * x; };

    for (const auto& [demoState, f_demo] : stateFrequencies) {
        double demoPayoff = 0.0;
        if (!payoffByState.empty()) {
            auto it = payoffByState.find(demoState);
            if (it != payoffByState.end()) demoPayoff = it->second;
        } else {
            for (size_t j = 0; j < demoState.size(); ++j) {
                if (demoState[j] == 1.0) demoPayoff += payoffs[j];
            }
        }

        const double prestigeBias = (slope == 0.0 ? 1.0 : std::pow(demoPayoff, slope));
        const double base = f_demo * prestigeBias;

        // Normalize within this demonstrator so the average expression multiplier is ~1
        double sumSq = 0.0;
        size_t k = 0;
        for (size_t t = 0; t < n; ++t) {
            if (demoState[t] == 0.0) continue;
            const double v = (dist == TraitDistribution::Depth) ? depths[t] : payoffs[t];
            sumSq += sq(v);
            ++k;
        }
        const double norm = (k > 0 && sumSq > 0.0) ? (sumSq / static_cast<double>(k)) : 1.0;

        for (size_t trait = 0; trait < n; ++trait) {
            if (repertoire[trait] != 0.0) continue; // learner already has this trait
            if (demoState[trait] != 1.0) continue; // demonstrator doesn't express the trait

            // Expression bias: scale by trait attribute and normalize within demonstrator
            const double v = (dist == TraitDistribution::Depth) ? depths[trait] : payoffs[trait];
            const double expr = sq(v) / norm;

            weights[trait] += base * expr;
        }
    }

    return weights;
}

std::vector<double> prestige2BaseWeights(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<double>& payoffs,
    const std::vector<double>& depths,
    TraitDistribution dist,
    double slope
) {
    // Prestige-by-size: sum over demonstrators who have trait i of f_demo * (|demo repertoire|^slope)
    std::vector<double> weights(repertoire.size(), 0.0);

    const size_t n = repertoire.size();

    // Default case: traits are expressed with uniform probability
    if (dist == TraitDistribution::Uniform) {
        for (size_t trait = 0; trait < n; ++trait) {
            if (repertoire[trait] != 0) continue; // learner already has this trait

            for (const auto& [demoState, f_demo] : stateFrequencies) {
                if (demoState[trait] != 1) continue; // demonstrator lacks the trait

                size_t demoSize = 0;
                for (double i : demoState) {
                    if (i == 1) ++demoSize;
                }

                weights[trait] += f_demo * std::pow(static_cast<double>(demoSize), slope);
            }
        }

        return weights;
    }

    auto sq = [](double x) { return x * x; };

    for (const auto& [demoState, f_demo] : stateFrequencies) {
        size_t demoSize = 0;
        for (double i : demoState) {
            if (i == 1) ++demoSize;
        }

        const double base = f_demo * std::pow(static_cast<double>(demoSize), slope);

        // Normalize within this demonstrator so the average expression multiplier is ~1
        double sumSq = 0.0;
        size_t k = 0;
        for (size_t t = 0; t < n; ++t) {
            if (demoState[t] == 0) continue;
            const double v = (dist == TraitDistribution::Depth) ? depths[t] : payoffs[t];
            sumSq += sq(v);
            ++k;
        }
        const double norm = (k > 0 && sumSq > 0.0) ? (sumSq / static_cast<double>(k)) : 1.0;

        for (size_t trait = 0; trait < n; ++trait) {
            if (repertoire[trait] != 0) continue; // learner already has this trait
            if (demoState[trait] != 1) continue; // demonstrator lacks the trait

            // Expression bias: scale by trait attribute and normalize within demonstrator
            const double v = (dist == TraitDistribution::Depth) ? depths[trait] : payoffs[trait];
            const double expr = sq(v) / norm;

            weights[trait] += base * expr;
        }
    }

    return weights;
}

std::vector<double> conformityBaseWeights(
    const std::vector<double>& traitFrequencies,
    double slope
) {
    std::vector<double> weights(traitFrequencies.size());
    
    // Conformity: frequency^slope 
    std::ranges::transform(traitFrequencies, weights.begin(), 
        [slope](double f) {
            return std::pow(f, slope);
        });
    
    return weights;
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
    const std::vector<double>& learnableProbs,
    const std::vector<double>& depths,
    TraitDistribution traitDist   
) {
    switch (strategy) {
    case Random:
        return traitFrequencies;
    case Payoff:
        return payoffBaseWeights(payoffs, traitFrequencies, slope);
    case Proximal:
        return proximalBaseWeights(repertoire, stateFrequencies, allStates, payoffs, depths, traitDist, slope);
    case Prestige:
        return prestigeBaseWeights(repertoire, stateFrequencies, allStates, statePayoffs, payoffs, depths, traitDist, slope);
    case Conformity:
        return conformityBaseWeights(traitFrequencies, slope);
    case Prestige2:
        return prestige2BaseWeights(repertoire, stateFrequencies,  payoffs, depths, traitDist, slope);
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
    double transparency,
    const std::vector<double>& statePayoffs,
    const std::vector<double>& learnableProbs,
    const AdjacencyMatrix& adjMatrix,
    const std::vector<double>& depths,
    TraitDistribution dist
)  {
    // Trait weights (typically w_i = f_i * bias_i)
    std::vector<double> weights = baseWeights(
        strategy, repertoire, payoffs, traitFrequencies, 
        stateFrequencies, allStates, slope, statePayoffs, learnableProbs, depths, dist);

    DEBUG_PRINT(2, "Current repertoire:");
    if (DEBUG_LEVEL >= 2) {
        for (double i: repertoire) {
            std::cout << (i == 1.0 ? "1" : "0");
        }
        std::cout << '\n';
    }

    DEBUG_PRINT(2, "Base weights:");
    if (DEBUG_LEVEL >= 2) {
        for (size_t i = 0; i < weights.size(); ++i) {
            std::cout << "Trait " << i <<": " << weights[i] << '\n';
        }
    }

    // Candidate weights: only traits not in repertoie and observable
    std::vector<double> candidateWeights(repertoire.size());
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        candidateWeights[trait] = (repertoire[trait] == 1.0 || traitFrequencies[trait] == 0.0) ? 0.0 : weights[trait];
    }

    DEBUG_PRINT (2, "Weights after filtering learned/unobserved traits:");
    if (DEBUG_LEVEL >= 2) {
        for (size_t i = 0; i < candidateWeights.size(); ++i) {
            std::cout << "Trait " << i <<": " << candidateWeights[i] << '\n';
        }
    }

     // Perceived learnability: penalize traits with more missing prerequisites (controlled by transparency).
    // missingPrereqFrac is normalized to [0,1] so (1 - missingPrereqFrac) is a “closeness” score.
    auto missingPrereqFrac = normalizedMissingPrereqs(repertoire, adjMatrix, transparency);

    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        candidateWeights[trait] = candidateWeights[trait] * std::pow((1 - missingPrereqFrac[trait]), transparency);
    }
    
    const double total = std::accumulate(candidateWeights.begin(), candidateWeights.end(), 0.0);

    if (total == 0.0) {
        DEBUG_PRINT(1, "All weights are zero, returning zero vector for repertoire");
        return std::vector<double>(repertoire.size(), 0.0);
    }

    // Normalize by total weight
    std::ranges::transform(candidateWeights, candidateWeights.begin(),
        [total](double w) { return w / total; });

    DEBUG_PRINT(2, "Normalized weights:");
    if (DEBUG_LEVEL >= 2) {
        for (size_t i = 0; i < candidateWeights.size(); ++i) {
            std::cout << "Trait " << i <<": " << candidateWeights[i] << '\n';
        }
    };

    return candidateWeights;
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
    double transparency,
    const AdjacencyMatrix& adjMatrix,
    const std::vector<double>& statePayoffs,
    const std::vector<double>& depths,
    TraitDistribution traitDist
) {
    std::vector<Repertoire> newStates = retrieveBetterRepertoires(allStates, repertoire);
    std::vector<double> learnableProbs = learnability(repertoire, adjMatrix);
    std::vector<double> w = normalizedWeights(strategy, repertoire, payoffs,
        traitFrequencies, stateFrequencies, newStates, slope, transparency,
        statePayoffs, learnableProbs, adjMatrix, depths, traitDist);
    

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
    double transparency,
    const std::vector<double>& depths,
    TraitDistribution dist
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

            auto transitions = transitionFromState(strategy, r, payoffs,
                traitFrequencies, stateFrequencies, allStates, slope,
                transparency, adjMatrix, statePayoffs, depths, dist);
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