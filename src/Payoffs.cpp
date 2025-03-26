#include "Payoffs.hpp"
#include <algorithm>
#include <stdexcept>
#include <fstream>

PayoffVector generatePayoffs(const std::vector<int>& distances, double alpha, const std::vector<size_t>& shuffleSequence, int payoffDist = 0) {
    size_t n = distances.size();
    size_t non_root_count = n - 1;
    PayoffVector payoffs(n, 0.0);

    if (n == 121) {
        std::string filePath = "../data/payoffs_121.csv";
        std::ifstream file(filePath);
        
        if (!file.is_open()) {
            throw std::runtime_error("Could not open payoffs file: " + filePath);
        }
        
        std::string line;
        size_t index = 0;
        
        while (std::getline(file, line) && index < n) {
            try {
                payoffs[index++] = std::stod(line);
            } catch (const std::exception& e) {
                throw std::runtime_error("Error parsing payoff value at line " + std::to_string(index) + ": " + e.what());
            }
        }
        
        if (index < n) {
            throw std::runtime_error("Not enough payoff values in file: expected " + std::to_string(n) + ", got " + std::to_string(index));
        }
        
        return payoffs;
    }
    
    // Root trait always has zero payoff
    payoffs[0] = 0.0;
    
    if (alpha <= 0.0) {
        // Use payoffDist approach when alpha is 0
        std::vector<double> non_root_payoffs(non_root_count);
        
        if (payoffDist == 0) {
            // Original equal spacing approach
            double spacing = 2.0 / (non_root_count + 1);
            for (size_t i = 0; i < non_root_count; ++i) {
                non_root_payoffs[i] = spacing * (i + 1);
            }
        } else {
            // High/low value approach
            size_t highValueCount = std::min(static_cast<size_t>(payoffDist), non_root_count);
            size_t lowValueCount = non_root_count - highValueCount;
            
            // Set values to maintain mean of 1.0
            double lowValue = 0.2;
            double highValue;
            
            if (lowValueCount == 0) {
                highValue = 1.0;
            } else {
                highValue = (non_root_count - (lowValueCount * lowValue)) / highValueCount;
            }
            
            for (size_t i = 0; i < highValueCount; ++i) {
                non_root_payoffs[i] = highValue;
            }
            for (size_t i = highValueCount; i < non_root_count; ++i) {
                non_root_payoffs[i] = lowValue;
            }
        }
        
        // Apply shuffle sequence
        std::vector<double> shuffled_non_root_payoffs(non_root_count);
        for (size_t i = 0; i < non_root_count; ++i) {
            shuffled_non_root_payoffs[i] = non_root_payoffs[shuffleSequence[i]];
        }
        
        // Assign to payoffs vector
        size_t idx = 0;
        for (size_t trait = 1; trait < n; ++trait) {
            payoffs[trait] = shuffled_non_root_payoffs[idx++];
        }
    } else {
        // Order traits by distance (depth) when alpha > 0
        std::vector<size_t> traitIndices;
        for (size_t i = 1; i < n; ++i) {
            traitIndices.push_back(i);
        }
        
        // Sort traits by distance from root (ascending)
        std::ranges::sort(traitIndices,
            [&distances](size_t a, size_t b) {
                return distances[a] < distances[b];
            });
        
        // Assign equidistant payoffs based on depth ordering
        double spacing = 2.0 / (non_root_count + 1);
        for (size_t i = 0; i < traitIndices.size(); ++i) {
            payoffs[traitIndices[i]] = spacing * (i + 1);
        }
    }
    
    // Adjust non-root payoffs to maintain a mean of 1
    double sum = 0.0;
    size_t count = 0;
    for (size_t trait = 1; trait < n; ++trait) {
        sum += payoffs[trait];
        count++;
    }
    
    double mean = (count > 0) ? sum / count : 0.0;
    
    if (mean > 0.0) {
        // Scale all non-root payoffs to get a mean of 1.0
        double scaleFactor = 1.0 / mean;
        for (size_t trait = 1; trait < n; ++trait) {
            payoffs[trait] *= scaleFactor;
        }
    }
    
    return payoffs;
}