#this compiles the program on snellius
module load 2023
module load CMake/3.26.3-GCCcore-12.3.0
module load zlib/1.3.1  
mkdir build && cd build
cmake ..
make 