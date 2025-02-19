#include "Learning.hpp"
#include "Debug.hpp"
#include "Types.hpp"
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include <unordered_set>

#ifndef PAYOFF_SENSITIVITY
#define PAYOFF_SENSITIVITY 1.0
#endif

std::vector<bool> learnability(const Repertoire& repertoire, const Parents& parents) {
    std::vector<bool> learnable(repertoire.size());

    for (size_t trait = 0; trait < repertoire.size(); ++trait) {
        auto traitParents = parents[trait];

        // true if all parents are in the repertoire
        bool parent_product = std::accumulate(traitParents.begin(), traitParents.end(), true,
            [&repertoire](bool acc, Trait parent) {
                return acc && repertoire[parent];
            });

        // true if the trait is not in the repertoire and all parents are in the repertoire
        bool is_learnable = !repertoire[trait] && parent_product;

        learnable[trait] = is_learnable;
    }
    return learnable;
}



double computeDelta(const Repertoire& r, const Repertoire& s) {
    // count the number of traits that are present in target state s but not in current state r
    return std::inner_product(s.begin(), s.end(), r.begin(), 0.0,
        std::plus<>(), [](bool s_i, bool r_i) { return (s_i && !r_i) ? 1 : 0; });
}

std::vector<double> proximalBaseWeights(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope
) {
    std::vector<double> w_star(repertoire.size(), 0.0);
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        if (!repertoire[trait]) {
            for (const auto& state : allStates) {
                if (state[trait]) {
                    auto freqIt = stateFrequencies.find(state);
                    if (freqIt != stateFrequencies.end()) {
                        double delta = computeDelta(repertoire, state);
                        if (delta > 0) {
                            w_star[trait] += freqIt->second * std::pow(slope, -delta);;
                        }
                    }
                }
            }
        }
    }
    return w_star;
}

std::vector<double> prestigeBaseWeights(
    const Repertoire& repertoire,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    double slope
) {
    std::vector<double> w_star(repertoire.size(), 0.0);
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        if (!repertoire[trait]) {
            for (const auto& state : allStates) {
                if (state[trait]) {
                    auto freqIt = stateFrequencies.find(state);
                    if (freqIt != stateFrequencies.end()) {
                        double delta = computeDelta(repertoire, state);
                        if (delta > 0) {
                            w_star[trait] += freqIt->second * std::pow(slope, delta);;
                        }
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
    
    std::transform(traitFrequencies.begin(), traitFrequencies.end(), w_star.begin(), [slope](double f) {return std::pow(f, slope);});

    return w_star;
}

std::vector<double> perfectBaseWeights(
    const Repertoire& repertoire,
    const PayoffVector& payoffs,
    const Parents& parents
) {
    // always pick the learnable trait with the highest payoff
    auto learnable = learnability(repertoire, parents);
    double highestPayoff = 0.0;
    double bestTraitidx = 0;
    std::vector<double> w_star(repertoire.size(), 0.0);
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        if (learnable[trait] && payoffs[trait] > highestPayoff) {
            highestPayoff = payoffs[trait];
            bestTraitidx = trait;
        }

    }
    w_star[bestTraitidx] = 1.0;
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
    const Parents& parents   
) {
    switch (strategy) {
    case Random:
        return traitFrequencies;
    case Payoff:
        {
            std::vector<double> result(payoffs.size());
            std::transform(payoffs.begin(), payoffs.end(), traitFrequencies.begin(), result.begin(),
               [slope](double payoff, double traitFrequency) {
                   return traitFrequency * std::pow(payoff, slope);
               });
            return result;
        }
    case Proximal:
        return proximalBaseWeights(repertoire, stateFrequencies, allStates, slope);
    case Prestige:
        return prestigeBaseWeights(repertoire, stateFrequencies, allStates, slope);
    case Conformity:
        return conformityBaseWeights(traitFrequencies, slope);
    case Perfect:
        return perfectBaseWeights(repertoire, payoffs, parents);
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
    const Parents& parents
)  {
    std::vector<double> w_star = baseWeights(strategy, repertoire, payoffs, traitFrequencies, stateFrequencies, allStates, slope, parents);

    DEBUG_PRINT(2, "Current repertoire:");
    if (DEBUG_LEVEL >= 2) {
        for (bool i : repertoire) {
            std::cout << (i ? "1" : "0");
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
        w_unlearned[trait] = (repertoire[trait] || traitFrequencies[trait] == 0.0) ? 0.0 : w_star[trait];
    }

    DEBUG_PRINT (2, "Weights after setting learned ones to 0:");
    if (DEBUG_LEVEL >= 2) {
        for (size_t i = 0; i < w_unlearned.size(); ++i) {
            std::cout << "Trait " << i <<": " << w_unlearned[i] << '\n';
        }
    }
    
    double total = std::accumulate(w_unlearned.begin(), w_unlearned.end(), 0.0);

    if (total == 0.0) {
        return std::vector<double>(repertoire.size(), 0.0);
    }

    std::transform(w_unlearned.begin(), w_unlearned.end(), w_unlearned.begin(),
        [total](double w) { return w / total; });
    


    DEBUG_PRINT(2, "Normalized weights:");
    if (DEBUG_LEVEL >= 2) {
        for (size_t i = 0; i < w_unlearned.size(); ++i) {
            std::cout << "Trait " << i <<": " << w_unlearned[i] << '\n';
        }
    };

    return w_unlearned;
}



Repertoire learnTrait(const Repertoire& repertoire, Trait trait) {
    Repertoire newRepertoire = repertoire;
    newRepertoire[trait] = true;
    return newRepertoire;
}

std::vector<std::pair<Repertoire, double>> transitionFromState(
    Strategy strategy,
    const Repertoire& repertoire, 
    const PayoffVector& payoffs, 
    const std::vector<double>& traitFrequencies,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,  
    const std::vector<Repertoire>& allStates,
    const Parents& parents,
    double slope
) {

    std::vector<Repertoire> newStates = retrieveBetterRepertoires(allStates, repertoire);
    std::vector<double> w = normalizedWeights(strategy, repertoire, payoffs, traitFrequencies, stateFrequencies, newStates, slope, parents);
    std::vector<bool> learnable = learnability(repertoire, parents);

    std::vector<std::pair<Repertoire, double>> transitions;

    for(Trait trait = 0; trait < repertoire.size(); ++trait) {
        if (learnable[trait] && w[trait] > 0.0) {
            transitions.emplace_back(learnTrait(repertoire, trait), w[trait]);
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

std::vector<Repertoire> generateReachableRepertoires(
    Strategy strategy, 
    const AdjacencyMatrix& adjMatrix, 
    const PayoffVector& payoffs, 
    const std::vector<double>& traitFrequencies,
    const std::unordered_map<Repertoire, double, RepertoireHash>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    const Parents& parents,
    double slope
) {
    size_t n = adjMatrix.size();
    Repertoire initialRepertoire(n, false);
    initialRepertoire[0] = true; // root trait is always learned

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

            auto transitions = transitionFromState(strategy, r, payoffs, traitFrequencies, stateFrequencies, allStates, parents, slope);
            for (const auto& transition : transitions) {
                const Repertoire& r_prime = transition.first;
                if (visited.find(r_prime) == visited.end()) {
                    queue.push(r_prime);
                }
            }
        }
    }

    return result;
}

std::vector<Repertoire> generateAllRepertoires(const AdjacencyMatrix& adjMatrix, const Parents& parents) {
    size_t n = adjMatrix.size();
    Repertoire initialRepertoire(n, false);
    initialRepertoire[0] = true; // root trait is always learned

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

            std::vector<bool> learnable = learnability(r, parents);
            for (Trait trait = 0; trait < n; ++trait) {
                if (learnable[trait]) {
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

size_t countLearnedTraits(const Repertoire& r) {
    return std::count(r.begin(), r.end(), true);
}

std::vector<Repertoire> retrieveBetterRepertoires(const std::vector<Repertoire>& repertoires, const Repertoire& singleRepertoire) {
    std::vector<Repertoire> result;
    for (const Repertoire& r : repertoires) {
        for (size_t trait = 0; trait < r.size(); ++trait) {
            if (r[trait] && !singleRepertoire[trait]) {
                result.push_back(r);
                break;
            }
        }
    }
    return result;
}