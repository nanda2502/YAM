#ifndef DEBUG_HPP
#define DEBUG_HPP

#include <iostream>


#define DEBUG_LEVEL 0  
std::ostream& setPrecision(std::ostream& os, int precision);

#define DEBUG_PRINT(level, x) \
    if (DEBUG_LEVEL >= level) { \
        std::ios_base::fmtflags old_flags = std::cout.flags(); \
        std::streamsize old_prec = std::cout.precision(); \
        std::cout << x << '\n'; \
        std::cout.flags(old_flags); \
        std::cout.precision(old_prec); \
    }

#endif //DEBUG_HPP
