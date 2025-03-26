#!/bin/bash
#SBATCH -p rome
#SBATCH -n 1
#SBATCH --cpus-per-task 5
#SBATCH -t 24:00:00

export OMP_NUM_THREADS=5

cd build

./yam 0 121



