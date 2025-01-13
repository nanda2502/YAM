#!/bin/bash
#SBATCH -p rome
#SBATCH -n 554
#SBATCH --cpus-per-task 4
#SBATCH -t 00:35:00

export OMP_NUM_THREADS=4

cd build

for i in {0..553}; do
    ./yam "$i" &
done

wait

