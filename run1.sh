#!/bin/bash
#SBATCH -p rome
#SBATCH -n 1
#SBATCH --cpus-per-task 128
#SBATCH -t 01:00:00

export OMP_NUM_THREADS=128

cd build

./yam 0 8&

wait
cd ..
./combine.sh



