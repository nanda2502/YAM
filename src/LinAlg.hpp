#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <vector>

std::pair<std::vector<std::vector<double>>, std::vector<int>> decomposeLU(const std::vector<std::vector<double>>& a);
std::vector<double> solveUsingLU(const std::vector<std::vector<double>>& LU, const std::vector<int>& P, const std::vector<double>& b);

#endif // LINEAR_ALGEBRA_HPP