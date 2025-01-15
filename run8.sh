#!/bin/bash
#SBATCH -p rome
#SBATCH -n 45
#SBATCH --cpus-per-task 32
#SBATCH -t 00:60:00

export OMP_NUM_THREADS=32

cd build

for i in {0..553}; do
    if [ ! -f "../output/expected_steps_${i}.csv.gz" ]; then
        ./yam "$i" 8&
    fi
done

wait



