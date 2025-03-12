#!/bin/bash
#SBATCH -p genoa
#SBATCH -n 1
#SBATCH --cpus-per-task 192
#SBATCH -t 01:00:00

export OMP_NUM_THREADS=96

cd build

./yam  50



