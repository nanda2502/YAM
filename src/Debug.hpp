#ifndef DEBUG_HPP
#define DEBUG_HPP

#define DEBUG_LEVEL 0  // Set to 0 to disable debug output, 1 to enable
#define DEBUG_PRINT(level, x) if (DEBUG_LEVEL >= level) { std::cout << x << '\n'; }

#endif //DEBUG_HPP
