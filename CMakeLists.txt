cmake_minimum_required(VERSION 3.26)
project(yam)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O2 -fopenmp")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "C:/msys64/mingw64")

include_directories(${CMAKE_SOURCE_DIR}/src)


# Main executable
file(GLOB_RECURSE MAIN_SOURCES
	"src/*.cpp"
	"src/markov/*.cpp"
	"src/learning/*.cpp"
)
add_executable(yam ${MAIN_SOURCES})




