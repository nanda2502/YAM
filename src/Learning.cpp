#include "Learning.hpp"
#include "Graph.hpp"
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <random>

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

std::vector<double> baseWeights(Strategy strategy, const PayoffVector& payoffs) {
    switch (strategy) {
    case RandomLearning:
    //vector of ones 
        return std::vector<double>(payoffs.size(), 1.0);
    case PayoffBasedLearning:
        return payoffs;
    default:
        throw std::runtime_error("Unknown strategy");
    }
}


std::vector<double> normalizedWeights(Strategy strategy, const Repertoire& repertoire, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies) {
    std::vector<double> w_star = baseWeights(strategy, payoffs);

    std::vector<double> w_unlearned(repertoire.size());
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        w_unlearned[trait] = repertoire[trait] ? 0.0 : w_star[trait] * traitFrequencies[trait];
    }

    double total = std::accumulate(w_unlearned.begin(), w_unlearned.end(), 0.0);

    if (total == 0.0) {
        return std::vector<double>(repertoire.size(), 0.0);
    }

    std::transform(w_unlearned.begin(), w_unlearned.end(), w_unlearned.begin(),
        [total](double w) { return w / total; });

    return w_unlearned;
}

std::vector<double> normalizedWeights(Strategy strategy, const Repertoire& repertoire, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies, std::mt19937& gen) {
    std::vector<double> w_star = baseWeights(strategy, payoffs);
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<double> w_unlearned(repertoire.size());
    for (Trait trait = 0; trait < repertoire.size(); ++trait) {
        if (repertoire[trait]) {
            w_unlearned[trait] = 0.0;
        } else {
            // Draw a random number and compare it to the trait frequency
            double random_value = dis(gen);
            w_unlearned[trait] = (random_value < traitFrequencies[trait]) ? w_star[trait] : 0.0;
        }
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

std::vector<std::pair<Repertoire, double>> transitionFromState(Strategy strategy, const Repertoire& repertoire, const AdjacencyMatrix& adjMatrix, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies) {
    //std::random_device rd;
    //std::mt19937 gen(rd());
    
    //std::vector<double> w = normalizedWeights(strategy, repertoire, payoffs, traitFrequencies, gen);
    std::vector<double> w = normalizedWeights(strategy, repertoire, payoffs, traitFrequencies);
    std::vector<bool> learnable = learnability(repertoire, adjMatrix);

    std::vector<std::pair<Repertoire, double>> transitions;

    for(Trait trait = 0; trait < repertoire.size(); ++trait) {
        if (learnable[trait] && w[trait] > 0.0) {
            transitions.emplace_back(learnTrait(repertoire, trait), w[trait]);
        }
    }

    return transitions;
}

double stayProbability(Strategy strategy, const Repertoire& repertoire, const AdjacencyMatrix& adjMatrix, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies) {
    auto transitions = transitionFromState(strategy, repertoire, adjMatrix, payoffs, traitFrequencies);

    double totalTransitionProbability = std::accumulate(
        transitions.begin(), transitions.end(), 0.0,
        [](double sum, const auto& transition) { 
            return sum + transition.second;
        }
    );

    return 1.0 - totalTransitionProbability;
}

std::vector<Repertoire> generateReachableRepertoires(Strategy strategy, const AdjacencyMatrix& adjMatrix, const PayoffVector& payoffs, const std::vector<double>& traitFrequencies) {
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

            auto transitions = transitionFromState(strategy, r, adjMatrix, payoffs, traitFrequencies);
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