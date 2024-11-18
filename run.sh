#!/bin/bash
#SBATCH -p rome
#SBATCH -n 150
#SBATCH --cpus-per-task=4
#SBATCH -t 00:15:00

export OMP_NUM_THREADS=4

cd build

for i in {0..149}; do
    ./yam "$i" &
done

wait