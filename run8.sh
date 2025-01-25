#!/bin/bash
#SBATCH -p rome
#SBATCH -n 554
#SBATCH --cpus-per-task 8
#SBATCH -t 00:30:00

export OMP_NUM_THREADS=8

cd build

for i in {0..553}; do
    if [ ! -f "../output/expected_steps_${i}.csv.gz" ]; then
        ./yam "$i" 8&
    fi
done

wait
cd ..
./combine.sh