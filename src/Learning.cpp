#include "Learning.hpp"
#include "Types.hpp"
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <queue>
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

double phi_proximal(double delta) {
    return std::pow(2.0, 1.0 - delta);
}

double phi_prestige(double delta) {
    return std::pow(2.0, delta - 1.0);
}

std::vector<double> proximalBaseWeights(
    const Repertoire& repertoire,
    const std::vector<double>& stateFrequencies,
    const std::vector<Repertoire>& allStates
) {
    std::vector<double> w_star(repertoire.size(), 0.0);
    // Compute trait contributions
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        if (!repertoire[trait]) {  // Only consider unlearned traits
            for (size_t stateidx = 0; stateidx < allStates.size(); ++stateidx) {
                const auto& state = allStates[stateidx];
                if (state[trait]) {  // If the state has trait j learned
                    double delta = computeDelta(repertoire, state);
                    if (delta > 0) {
                        int N_s = delta;  // Number of additional traits known in s but not in r
                        w_star[trait] += stateFrequencies[stateidx] * (phi_proximal(delta) / N_s); 
                    }
                }
            }
        }
    }
    return w_star;
}

std::vector<double> prestigeBaseWeights(
    const Repertoire& repertoire,
    const std::vector<double>& stateFrequencies,
    const std::vector<Repertoire>& allStates
) {
    std::vector<double> w_star(repertoire.size(), 0.0);
    // Compute trait contributions
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        if (!repertoire[trait]) {  // Only consider unlearned traits
            for (size_t stateidx = 0; stateidx < allStates.size(); ++stateidx) {
                const auto& state = allStates[stateidx];
                if (state[trait]) {  // If the state has trait j learned
                    double delta = computeDelta(repertoire, state);
                    if (delta > 0) {
                        int N_s = delta;  // Number of additional traits known in s but not in r
                        w_star[trait] += stateFrequencies[stateidx] * (phi_prestige(delta) / N_s); 
                    }
                }
            }
        }
    }
    return w_star;
}

double S_curve(double x, double total) {
    // between 0 and 1;
    return 1/(1+std::exp(-5 * ((x/total)-0.5)));
}


std::vector<double> conformityBaseWeights(
    const std::vector<double>& traitFrequencies,
    const Repertoire& repertoire
) {
    double total = 0.0;

    //Only count unlearned traits, so that the normalization works
    for (size_t i = 0; i < traitFrequencies.size(); i++) {
        total += repertoire[i] ? 0.0 : traitFrequencies[i];
    }

    std::vector<double> w_star(traitFrequencies.size());
    
    std::transform(traitFrequencies.begin(), traitFrequencies.end(), w_star.begin(), [total](double f) {return S_curve(f, total);});

    return w_star;
}

std::vector<double> baseWeights(
    Strategy strategy,
    const Repertoire& repertoire,
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::vector<double>& stateFrequencies,
    const std::vector<Repertoire>& allStates
) {
    switch (strategy) {
    case RandomLearning:
        return traitFrequencies;
    case PayoffBasedLearning:
        {
            std::vector<double> result(payoffs.size());
            std::transform(payoffs.begin(), payoffs.end(), traitFrequencies.begin(), result.begin(),
               [](double payoff, double traitFrequency) {
                   return traitFrequency * std::pow(payoff, 5);
               });
            return result;
        }
    case ProximalLearning:
        return proximalBaseWeights(repertoire, stateFrequencies, allStates);
    case PrestigeBasedLearning:
        return prestigeBaseWeights(repertoire, stateFrequencies, allStates);
    case ConformityBasedLearning:
        return conformityBaseWeights(traitFrequencies, repertoire);
    default:
        throw std::runtime_error("Unknown strategy");
    }
}


std::vector<double> normalizedWeights(
    Strategy strategy,
    const Repertoire& repertoire,
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::vector<double>& stateFrequencies,
    const std::vector<Repertoire>& allStates
)  {
    std::vector<double> w_star = baseWeights(strategy, repertoire, payoffs, traitFrequencies, stateFrequencies, allStates);

    std::vector<double> w_unlearned(repertoire.size());
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        w_unlearned[trait] = repertoire[trait] ? 0.0 : w_star[trait];
    }

    double total = std::accumulate(w_unlearned.begin(), w_unlearned.end(), 0.0);

    if (total == 0.0) {
        return std::vector<double>(repertoire.size(), 0.0);
    }

    std::transform(w_unlearned.begin(), w_unlearned.end(), w_unlearned.begin(),
        [total](double w) { return w / total; });

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
    const std::vector<double>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    const Parents& parents
) {

    std::vector<Repertoire> newStates = retrieveBetterRepertoires(allStates, repertoire);
    std::vector<double> w = normalizedWeights(strategy, repertoire, payoffs, traitFrequencies, stateFrequencies, newStates);
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
    const std::vector<double>& stateFrequencies,
    const std::vector<Repertoire>& allStates,
    const Parents& parents
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

            auto transitions = transitionFromState(strategy, r, payoffs, traitFrequencies, stateFrequencies, allStates, parents);
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
    size_t singleRepertoireLearnedTraits = countLearnedTraits(singleRepertoire);

    for (const Repertoire& r : repertoires) {
        if (countLearnedTraits(r) > singleRepertoireLearnedTraits) {
            result.push_back(r);
        }
    }

    return result;
}