#!/bin/bash
#SBATCH -p genoa
#SBATCH -N 1
#SBATCH --cpus-per-task 192
#SBATCH -t 24:00:00

export OMP_NUM_THREADS=192

cd build

./yam 0 8



