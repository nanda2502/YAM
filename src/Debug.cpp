#include "Debug.hpp"

std::ostream& setPrecision(std::ostream& os, int precision) {
    os.setf(std::ios::fixed);
    os.precision(precision);
    return os;
}