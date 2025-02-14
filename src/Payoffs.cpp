#include "Payoffs.hpp"



PayoffVector generateRawPayoffs(const std::vector<int>& distances, const std::vector<size_t>& shuffleSequence) {
    constexpr double mean = 1.0;
    size_t n = distances.size();
    size_t non_root_count = n - 1; 

    double spacing = (2.0 * mean) / (non_root_count + 1);
    std::vector<double> non_root_payoffs(non_root_count);

    // Generate equally spaced payoffs for non-root traits
    for (size_t i = 0; i < non_root_count; ++i) {
        non_root_payoffs[i] = spacing * (i + 1);
    }

    // Create a new vector for shuffled payoffs
    std::vector<double> shuffled_non_root_payoffs(non_root_count);

    // Reorder the non_root_payoffs according to the shuffle sequence
    for (size_t i = 0; i < non_root_count; ++i) {
        shuffled_non_root_payoffs[i] = non_root_payoffs[shuffleSequence[i]];
    }

    PayoffVector payoffs(n, 0.0);

    size_t non_root_idx = 0;
    for (size_t trait = 1; trait < n; ++trait) {
        payoffs[trait] = shuffled_non_root_payoffs[non_root_idx++];
    }

    return payoffs;
}

void addDistanceBonuses(PayoffVector& payoffs, const std::vector<int>& distances, double alpha) {
    double original_mean = 0.0;
    size_t non_root_count = 0;
    
    // Calculate original mean of non-root payoffs
    for (size_t trait = 1; trait < payoffs.size(); ++trait) {
        original_mean += payoffs[trait];
        non_root_count++;
    }
    original_mean /= non_root_count;

    // Add distance bonuses
    for (size_t trait = 0; trait < payoffs.size(); ++trait) {
        if (distances[trait] != 0) {
            payoffs[trait] += alpha * distances[trait];
        }
    }

    // Calculate new mean after adding bonuses
    double new_mean = 0.0;
    for (size_t trait = 1; trait < payoffs.size(); ++trait) {
        new_mean += payoffs[trait];
    }
    new_mean /= non_root_count;

    // Scale payoffs to maintain original mean
    double scaling_factor = original_mean / new_mean;
    for (size_t trait = 0; trait < payoffs.size(); ++trait) {
        payoffs[trait] *= scaling_factor;
    }
}

PayoffVector generatePayoffs(const std::vector<int>& distances, double alpha, const std::vector<size_t>& shuffleSequence) {
    PayoffVector payoffs = generateRawPayoffs(distances, shuffleSequence);
    addDistanceBonuses(payoffs, distances, alpha);
    return payoffs;
}
