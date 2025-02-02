#include "Payoffs.hpp"

PayoffVector generateRawPayoffs(const std::vector<size_t>& shuffleSequence) {
    constexpr double mean = 1.0;
    constexpr size_t total_payoffs = 8; // Always generate 8 payoffs
    size_t non_root_count = total_payoffs - 1; // 7 non-root payoffs

    double spacing = (2.0 * mean) / (non_root_count + 1);
    std::vector<double> non_root_payoffs(non_root_count);

    // Generate equally spaced payoffs for non-root traits
    for (size_t i = 0; i < non_root_count; ++i) {
        non_root_payoffs[i] = spacing * (i + 1);
    }

    // Determine the number of payoffs to generate based on the shuffle sequence length
    size_t num_payoffs = shuffleSequence.size();

    // Create a new vector for shuffled payoffs
    std::vector<double> shuffled_non_root_payoffs(num_payoffs);

    // Reorder the non_root_payoffs according to the shuffle sequence
    for (size_t i = 0; i < num_payoffs; ++i) {
        // Ensure that the shuffle sequence index is within bounds
        if (shuffleSequence[i] < non_root_count) {
            shuffled_non_root_payoffs[i] = non_root_payoffs[shuffleSequence[i]];
        } else {
            // Handle out-of-bounds access (e.g., set to 0 or some default value)
            shuffled_non_root_payoffs[i] = 0.0;
        }
    }

    // Create the payoff vector with the same length as the shuffle sequence
    PayoffVector payoffs(num_payoffs + 1, 0.0); // +1 for the root payoff

    // Assign the shuffled payoffs to the non-root traits
    for (size_t i = 0; i < num_payoffs; ++i) {
        payoffs[i + 1] = shuffled_non_root_payoffs[i]; // +1 to skip the root
    }

    return payoffs;
}

void addDistanceBonuses(PayoffVector& payoffs, const std::vector<int>& distances, double alpha) {
    for (size_t trait = 0; trait < payoffs.size(); ++trait) {
        if (distances[trait] != 0) {
            payoffs[trait] += alpha * distances[trait];
        }
    }
}

PayoffVector generatePayoffs(const std::vector<int>& distances, double alpha, const std::vector<size_t>& shuffleSequence) {
    PayoffVector payoffs = generateRawPayoffs(shuffleSequence);
    addDistanceBonuses(payoffs, distances, alpha);
    return payoffs;
}
