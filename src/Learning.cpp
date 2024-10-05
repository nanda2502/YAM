#include "Learning.hpp"
#include "Graph.hpp"
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include <unordered_set>

std::vector<bool> learnability(const Repertoire& repertoire, const AdjacencyMatrix& adjMatrix) {
    std::vector<bool> learnable(repertoire.size());

    for (size_t trait = 0; trait < repertoire.size(); ++trait) {
        auto parents = parentTraits(adjMatrix, trait);

        // true if all parents are in the repertoire
        bool parent_product = std::accumulate(parents.begin(), parents.end(), true,
            [&repertoire](bool acc, Trait parent) {
                return acc && repertoire[parent];
            });

        // true if the trait is not in the repertoire and all parents are in the repertoire
        bool is_learnable = !repertoire[trait] && parent_product;

        learnable[trait] = is_learnable;
    }
    return learnable;
}

std::unordered_map<Repertoire, double, RepertoireHash> computeStateFrequencies(
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates) {

    std::unordered_map<Repertoire, double, RepertoireHash> stateFrequencies;

    for (const auto& state : allStates) {
        double probability = 1.0;
        for (size_t i = 0; i < state.size(); ++i) {
            if (state[i]) {
                probability *= traitFrequencies[i]; // Probability trait i is learned
            } else {
                probability *= (1 - traitFrequencies[i]); // Probability trait i is not learned
            }
        }
        stateFrequencies[state] = probability;
    }

    return stateFrequencies;
}


double computeDelta(const Repertoire& r, const Repertoire& s) {
    // count the number of traits that are present in target state s but not in current state r
    return std::inner_product(s.begin(), s.end(), r.begin(), 0.0,
        std::plus<>(), [](bool s_i, bool r_i) { return (s_i && !r_i) ? 1 : 0; });
}

double phi(double delta) {
    return std::pow(2.0, 1.0 - delta);
}

std::vector<double> vygotskianBaseWeights(
    const Repertoire& repertoire,
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates
) {
    std::vector<double> w_star(repertoire.size(), 0.0);
    auto stateFrequencies = computeStateFrequencies(traitFrequencies, allStates);
    // Compute trait contributions
    for (size_t j = 0; j < repertoire.size(); ++j) {
        if (!repertoire[j]) {  // Only consider unlearned traits
            for (const auto& state : allStates) {
                if (state[j]) {  // If the state has trait j learned
                    double delta = computeDelta(repertoire, state);
                    if (delta > 0) {
                        int N_s = delta;  // Number of additional traits known in s but not in r
                        w_star[j] += stateFrequencies.at(state) / N_s;
                    }
                }
            }
        }
    }

    return w_star;
}

std::vector<double> baseWeights(
    Strategy strategy,
    const Repertoire& repertoire,
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates  
) {
    switch (strategy) {
    case RandomLearning:
        return traitFrequencies;
    case PayoffBasedLearning:
        {
            std::vector<double> result(payoffs.size());
            std::transform(payoffs.begin(), payoffs.end(), traitFrequencies.begin(), result.begin(), std::multiplies<>());
            return result;
        }
    case VygotskianLearning:
        return vygotskianBaseWeights(repertoire, traitFrequencies, allStates);
    default:
        throw std::runtime_error("Unknown strategy");
    }
}


std::vector<double> normalizedWeights(
    Strategy strategy,
    const Repertoire& repertoire,
    const PayoffVector& payoffs,
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates
)  {
    std::vector<double> w_star = baseWeights(strategy, repertoire, payoffs, traitFrequencies, allStates);

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
    const AdjacencyMatrix& adjMatrix, 
    const PayoffVector& payoffs, 
    const std::vector<double>& traitFrequencies,
    const std::vector<Repertoire>& allStates
) {

    std::vector<Repertoire> newStates = retrieveBetterRepertoires(allStates, repertoire);
    std::vector<double> w = normalizedWeights(strategy, repertoire, payoffs, traitFrequencies, newStates);
    std::vector<bool> learnable = learnability(repertoire, adjMatrix);

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
    const std::vector<Repertoire>& allStates
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

            auto transitions = transitionFromState(strategy, r, adjMatrix, payoffs, traitFrequencies, allStates);
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

std::vector<Repertoire> generateAllRepertoires(const AdjacencyMatrix& adjMatrix) {
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

            std::vector<bool> learnable = learnability(r, adjMatrix);
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