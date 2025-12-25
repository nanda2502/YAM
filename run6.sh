#!/bin/bash
#SBATCH -p rome
#SBATCH -n 57
#SBATCH --cpus-per-task 4
#SBATCH -t 00:10:00

export OMP_NUM_THREADS=4

cd build

for i in {0..56}; do
    if [ ! -f "../output/yam_out_${i}.csv" ]; then
        ./yam "$i" 6&
    fi
done

wait