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

std::vector<double> proximalBaseWeightsOld(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope
) {
    std::vector<double> w_star(repertoire.size(), 0.0);
    
    // Loop over each trait
    for (size_t trait = 0; trait < repertoire.size(); ++trait) {
        if (repertoire[trait] == 0) { // Agent doesn't have this trait
            // Loop over all potential demonstrator states in the stateFrequencies map
            for (const auto& [state, frequency] : stateFrequencies) {
                if (state[trait] == 1) { // Demonstrator has this trait
                    auto delta = computeDelta(repertoire, state);
                    if (delta > 0) {
                        // Add to the weight using inverse of delta
                        w_star[trait] += frequency * std::pow(delta, -slope);
                    }
                }
            }
        }
    }
    
    return w_star;
}

std::vector<double> proximalBaseWeights_meandist(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope
) {
    std::vector<double> w_star(repertoire.size(), 0.0);

    // learner's repertoire size
    int learnerSize = 0;
    for (double val : repertoire) {
        if (val == 1.0) learnerSize++;
    }

    // Precompute trait × size frequencies
    size_t numTraits = repertoire.size();
    size_t maxSize = numTraits;
    std::vector<std::vector<double>> traitSizeFreq(numTraits, std::vector<double>(maxSize + 1, 0.0));

    for (const auto& [state, frequency] : stateFrequencies) {
        int demoSize = 0;
        for (double val : state) {
            if (val == 1.0) demoSize++;
        }
        for (size_t trait = 0; trait < numTraits; ++trait) {
            if (state[trait] == 1.0) {
                traitSizeFreq[trait][demoSize] += frequency;
            }
        }
    }

    // Compute weights for each unlearned trait using size-based delta
    for (size_t trait = 0; trait < numTraits; ++trait) {
        if (repertoire[trait] == 0.0) {
            double sum = 0.0;
            for (int demoSize = learnerSize + 1; demoSize <= (int)maxSize; ++demoSize) {
                int delta = demoSize - learnerSize;
                if (delta > 0) {
                    sum += traitSizeFreq[trait][demoSize] * std::pow(delta, -slope);
                }
            }
            w_star[trait] = sum;
        }
    }

    return w_star;
}


std::vector<double> proximalBaseWeights(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope
) {
    const size_t T = repertoire.size();
    std::vector<double> w_star(T, 0.0);

    // Learner size m := |S(i)| = number of present traits in the focal repertoire
    int m = 0;
    for (double x : repertoire) if (x == 1.0) ++m;

    // Build mean-field summaries from the population measure 'stateFrequencies'.
    //
    // θ_u(k)  := Pr(u ∧ |S|=k) — population mass of demonstrators of size k carrying trait u.
    std::vector<std::vector<double>> theta(T, std::vector<double>(T + 1, 0.0));

    // ψ_num[u][v][k] accumulates Pr(u ∧ v ∧ |S|=k).  Later we divide by θ_u(k) to obtain
    // ψ_{u→v}(k) = Pr(v | u, |S|=k), i.e., the typical "bundle" of other traits that co-occur with u at size k.
    std::vector<std::vector<std::vector<double>>> psi_num(
        T, std::vector<std::vector<double>>(T, std::vector<double>(T + 1, 0.0))
    );

    // For the learner’s size m we also need g_v(m) = Pr(v | |S|=m),
    // the prevalence of each trait v among peers of the same size.
    std::vector<double> g_num(T, 0.0);
    double Zm = 0.0; // Zm = Pr(|S|=m) (normalizer for the size-m slice)

    // Single pass to fill θ, ψ numerators, and the size-m slice for g(·|m).
    for (const auto& kv : stateFrequencies) {
        const Repertoire& state = kv.first;
        const double freq = kv.second;

        int k = 0;
        for (double x : state) if (x == 1.0) ++k;

        // Contribute to g(·|m): restrict to the size-m subpopulation
        if (k == m) {
            Zm += freq;
            for (size_t v = 0; v < T; ++v) if (state[v] == 1.0) g_num[v] += freq;
        }

        // Contribute to θ and ψ numerators across all traits present in this state
        for (size_t u = 0; u < T; ++u) if (state[u] == 1.0) {
            theta[u][k] += freq;                       // adds to Pr(u ∧ |S|=k)
            for (size_t v = 0; v < T; ++v) if (state[v] == 1.0) {
                psi_num[u][v][k] += freq;              // adds to Pr(u ∧ v ∧ |S|=k)
            }
        }
    }

    // Convert g_num → g(·|m) = Pr(v | |S|=m)
    std::vector<double> g(T, 0.0);
    if (Zm > 0.0) {
        for (size_t v = 0; v < T; ++v) g[v] = g_num[v] / Zm;
    }

    // Small floor to keep the kernel Δ^{−slope} well-defined
    const double eps = 1e-12;

    // For each missing trait u, compute the proximal mean-field weight:
    //
    // w*(u) = Σ_k θ_u(k) · ϕ( E[Δ | m,u,k] ),   with ϕ(x)=x^{−slope}
    //
    // where the "expected directed gap" is
    //   E[Δ | m,u,k] = k − Σ_{v≠u} Pr(v | u, k) · Pr(v | m)
    //
    // Interpretation:
    //   • A size-k demonstrator with u carries ~k−1 other traits.
    //   • A size-m learner typically shares g(v|m) of each other trait v.
    //   • Subtract expected overlap to obtain the expected shortfall Δ against such demonstrators.
    for (size_t u = 0; u < T; ++u) {
        if (repertoire[u] != 0.0) continue; // evaluate only u ∉ S(i)

        double sum = 0.0;

        // Aggregate contributions from demonstrator size classes k
        for (size_t k = 0; k <= T; ++k) {
            const double mass_uk = theta[u][k]; // θ_u(k) = Pr(u ∧ |S|=k)
            if (mass_uk <= 0.0) continue;

            // Expected overlap on "other" traits:
            // shared_other = Σ_{v≠u} ψ_{u→v}(k) · g(v|m)
            double shared_other = 0.0;
            for (size_t v = 0; v < T; ++v) {
                if (v == u) continue;
                const double psi_cond = psi_num[u][v][k] / mass_uk; // Pr(v | u, k)
                shared_other += psi_cond * g[v];
            }

            // E[Δ | m,u,k] = k − shared_other
            const double expected_delta = static_cast<double>(k) - shared_other;

            // θ_u(k) · (E[Δ])^{−slope}
            if (expected_delta > eps) {
                sum += mass_uk * std::pow(expected_delta, -slope);
            }
        }

        w_star[u] = sum;
    }

    return w_star;
}

std::vector<double> payoffBaseWeights(
    const std::vector<double>& payoffs, 
    const std::vector<double>& traitFrequencies, 
    double slope
) {
    std::vector<double> w_star(payoffs.size());
    
    for (size_t i = 0; i < payoffs.size(); ++i) {
        w_star[i] = traitFrequencies[i] * std::pow(payoffs[i], slope);
    }
    
    return w_star;
}

std::vector<double> prestigeBaseWeights(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    const std::vector<double>& statePayoffs,
    const std::vector<double>& payoffs,
    double slope
) {
    std::vector<double> w_star(repertoire.size(), 0.0);
    
    // Loop over each trait
    for (size_t trait = 0; trait < repertoire.size(); ++trait) {
        if (repertoire[trait] == 0) { // Agent doesn't have this trait
            // Loop over all potential demonstrator states in the stateFrequencies map
            for (const auto& [state, frequency] : stateFrequencies) {
                if (state[trait] == 1) { // Demonstrator has this trait
                    auto delta = computeDelta(repertoire, state);
                    if (delta > 0) {
                        // Find the state in allStates to get its payoff
                        auto stateIt = std::find(allStates.begin(), allStates.end(), state);
                        if (stateIt != allStates.end()) {
                            size_t stateIndex = std::distance(allStates.begin(), stateIt);
                            //double statePayoff = statePayoffs[stateIndex];
                            double statePayoff = 0.0;
                            for (size_t j = 0; j < state.size(); ++j) {
                                if (state[j] == 1.0) {
                                    statePayoff += payoffs[j];
                                }
                            }
                            // Add to the weight using payoff bias
                            w_star[trait] += frequency * std::pow(statePayoff, slope);
                        }
                    }
                }
            }
        }
    }
    
    return w_star;
}

std::vector<double> prestige2BaseWeights(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope
) {
    // Used for Figure S3. This implementation weights by repertoire size instead of repertoire payoff
    std::vector<double> w_star(repertoire.size(), 0.0);
    
    // Loop over each trait
    for (size_t trait = 0; trait < repertoire.size(); ++trait) {
        if (repertoire[trait] == 0) { // Agent doesn't have this trait
            // Loop over all potential demonstrator states in the stateFrequencies map
            for (const auto& [state, frequency] : stateFrequencies) {
                if (state[trait] == 1) { // Demonstrator has this trait
                    auto delta = computeDelta(repertoire, state);
                    if (delta > 0) {
                        // Calculate repertoire size (number of traits the demonstrator has)
                        size_t repertoireSize = 0;
                        for (size_t i = 0; i < state.size(); ++i) {
                            if (state[i] == 1) {
                                repertoireSize++;
                            }
                        }
                        // Add to the weight using repertoire size bias
                        w_star[trait] += frequency * std::pow(repertoireSize, slope);
                    }
                }
            }
        }
    }
    
    return w_star;
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
    // This is not actually optimal because it doesn't consider how traits enable the learning of downstream traits.
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
        return prestigeBaseWeights(repertoire, stateFrequencies, allStates, statePayoffs, payoffs, slope);
    case Conformity:
        return conformityBaseWeights(traitFrequencies, slope);
    case Anticonformity:
        return anticonformityBaseWeights(traitFrequencies, slope);
    case Perfect:
        return perfectBaseWeights(repertoire, payoffs, learnableProbs);
    case Prestige2:
        return prestige2BaseWeights(repertoire, stateFrequencies, allStates, slope);
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

    // Attempt weights depend on partial structure knowledge. When transparency is 0, structure is completely opaque. 
    auto missingPrereqs = normalizedMissingPrereqs(repertoire, adjMatrix, transparency);

    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        w_unlearned[trait] = w_unlearned[trait] * std::pow((1 - missingPrereqs[trait]), transparency);
    }
    
    double total = std::accumulate(w_unlearned.begin(), w_unlearned.end(), 0.0);

    if (total == 0.0) {
        DEBUG_PRINT(1, "All weights are zero, returning zero vector for repertoire");
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
    double transparency,
    const AdjacencyMatrix& adjMatrix,
    const std::vector<double>& statePayoffs
) {
    std::vector<Repertoire> newStates = retrieveBetterRepertoires(allStates, repertoire);
    std::vector<double> learnableProbs = learnability(repertoire, adjMatrix);
    std::vector<double> w = normalizedWeights(strategy, repertoire, payoffs, traitFrequencies, stateFrequencies, newStates, slope, transparency, statePayoffs, learnableProbs, adjMatrix);
    

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
    double transparency
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

            auto transitions = transitionFromState(strategy, r, payoffs, traitFrequencies, stateFrequencies, allStates, slope, transparency, adjMatrix, statePayoffs);
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