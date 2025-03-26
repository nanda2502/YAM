#include "Learning.hpp"
#include "Debug.hpp"
#include "Types.hpp"
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <queue>
#include <unordered_map>
#include <unordered_set>



std::vector<double> learnability(
    const Repertoire& repertoire,
    const Parents& parents,
    double edgeWeight //represents the probability that a node with an edge to another node is a parent node of that node
) {
    std::vector<double> learnable(repertoire.size());

    for (size_t trait = 0; trait < repertoire.size(); ++trait) {
        auto traitParents = parents[trait];

        // Calculate probability that all parents are in the repertoire
        double parent_product = std::accumulate(traitParents.begin(), traitParents.end(), 1.0,
            [&repertoire, edgeWeight](double acc, Trait parent) {
                // If parent is in repertoire, multiply by edgeWeight (probability it's actually a parent)
                // If parent is not in repertoire, multiply by (1 - edgeWeight) (probability it's not a parent)
                return acc * (repertoire[parent] == 1.0 ? edgeWeight : (1.0 - edgeWeight));
            });

        // Trait is learnable if it's not in the repertoire and all parents are in the repertoire
        learnable[trait] = (repertoire[trait] == 0.0) ? parent_product : 0.0;
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
    std::vector<double> w_star(repertoire.size(), 0.0);
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        if (repertoire[trait] == 0.0) {
            for (const auto& state : allStates) {
                if (state[trait] == 1.0) {
                    auto it = stateFrequencies.find(state);
                    if (it != stateFrequencies.end()) {
                        auto delta = computeDelta(repertoire, state);
                        if (delta > 0) {
                            w_star[trait] += it->second * std::pow(delta, -slope);
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
    
    // First loop over all potential demonstrator states
    for (const auto& state : allStates) {
        auto it = stateFrequencies.find(state);
        if (it != stateFrequencies.end()) {
            auto delta = computeDelta(repertoire, state);
            if (delta > 0) {  // State has at least one trait not in repertoire
                // Count total traits in the demonstrator state
                int totalTraits = std::count(state.begin(), state.end(), 1.0);
                // Weight using total traits
                double stateWeight = it->second * std::pow(totalTraits, slope);
                
                // Then loop over traits to assign weights
                for (Trait trait = 0; trait < repertoire.size(); ++trait) {
                    if (repertoire[trait] == 0.0 && state[trait] == 1.0) {  // Trait is unlearned by agent but present in demonstrator
                        w_star[trait] += stateWeight;
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
    
    std::ranges::transform(traitFrequencies, w_star.begin(), [slope](double f) {return std::pow(f, slope);});

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

std::vector<double> payoffBaseWeights(const std::vector<double>& payoffs, const std::vector<double>& traitFrequencies, double slope) {
    std::vector<double> w_star(payoffs.size());
    std::ranges::transform(payoffs, traitFrequencies, w_star.begin(),
        [slope](double payoff, double traitFrequency) {
            return traitFrequency * std::pow(payoff, slope);
        });
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
        return payoffBaseWeights(payoffs, traitFrequencies, slope);
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

    std::ranges::transform(w_unlearned, w_unlearned.begin(),
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
    const Parents& parents,
    double slope,
    double edgeWeight
) {
    std::vector<Repertoire> newStates = retrieveBetterRepertoires(allStates, repertoire);
    std::vector<double> w = normalizedWeights(strategy, repertoire, payoffs, traitFrequencies, stateFrequencies, newStates, slope, parents);
    std::vector<double> learnableProbs = learnability(repertoire, parents, edgeWeight);

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
    const Parents& parents,
    double slope
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

            auto transitions = transitionFromState(strategy, r, payoffs, traitFrequencies, stateFrequencies, allStates, parents, slope);
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

std::vector<Repertoire> generateAllRepertoires(const AdjacencyMatrix& adjMatrix, const Parents& parents) {
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