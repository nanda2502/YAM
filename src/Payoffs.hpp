#ifndef PAYOFFS_HPP
#define PAYOFFS_HPP

#include "Types.hpp"
#include <vector>

PayoffVector generatePayoffs(const std::vector<double>& distances, double alpha, const std::vector<size_t>& shuffleSequence, int payoffDist);

#endif // PAYOFFS_HPP