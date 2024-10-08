#include "LinAlg.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

const double EPSILON = 1e-10;

std::pair<std::vector<std::vector<double>>, std::vector<int>> decomposeLU(const std::vector<std::vector<double>>& a) {
    int n = a.size();
    std::vector<std::vector<double>> LU = a;
    std::vector<int> P(n);
    for (int i = 0; i < n; i++) P[i] = i;

    for (int k = 0; k < n - 1; k++) {
        // Find the pivot row
        int pivot_row = k;
        double pivot_val = std::abs(LU[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (std::abs(LU[i][k]) > pivot_val) {
                pivot_row = i;
                pivot_val = std::abs(LU[i][k]);
            }
        }

        if (pivot_val < EPSILON) {
            throw std::runtime_error("Matrix is singular or nearly singular");
        }

        if (pivot_row != k) {
            std::swap(LU[k], LU[pivot_row]);
            std::swap(P[k], P[pivot_row]);
        }

        for (int i = k + 1; i < n; i++) {
            LU[i][k] /= LU[k][k];
            for (int j = k + 1; j < n; j++) {
                LU[i][j] -= LU[i][k] * LU[k][j];
            }
        }
    }

    return {LU, P};
}

std::vector<double> solveUsingLU(const std::vector<std::vector<double>>& LU, const std::vector<int>& P, const std::vector<double>& b) {
    int n = LU.size();
    
    // Forward substitution to solve Ly = Pb
    std::vector<double> y(n);
    for (int i = 0; i < n; i++) {
        y[i] = b[P[i]];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i][j] * y[j];
        }
    }

    // Backward substitution to solve Ux = y
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= LU[i][j] * x[j];
        }
        x[i] /= LU[i][i];
    }

    return x;
}