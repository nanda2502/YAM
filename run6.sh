#!/bin/bash
#SBATCH -p rome
#SBATCH -n 56
#SBATCH --cpus-per-task 4
#SBATCH -t 00:05:00

export OMP_NUM_THREADS=4

cd build

for i in {0..56}; do
    if [ ! -f "../output/expected_steps_${i}.csv.gz" ]; then
        ./yam "$i" 6&
    fi
done

wait