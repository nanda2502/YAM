#ifndef PAYOFFS_HPP
#define PAYOFFS_HPP

#include "Types.hpp"
#include <vector>
#include <random>

PayoffVector generatePayoffs(const std::vector<int>& distances, double alpha, std::mt19937& gen);

#endif // PAYOFFS_HPP