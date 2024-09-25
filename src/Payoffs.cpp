#include "Payoffs.hpp"
#include <random>


PayoffVector generatePayoffs(const std::vector<int>& distances, double alpha, std::mt19937& gen) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    PayoffVector payoffs(distances.size());

    for (size_t i = 0; i < distances.size(); ++i) {
        if (distances[i] == 0) {
            payoffs[i] = 0.0;
        } else {
            double randValue = dist(gen);
            double payoff = randValue + alpha * distances[i];
            payoffs[i] = payoff;
        }
    }

    return payoffs;
}