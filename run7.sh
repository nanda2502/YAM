#!/bin/bash
#SBATCH -p rome
#SBATCH -n 150
#SBATCH --cpus-per-task 4
#SBATCH -t 00:30:00

export OMP_NUM_THREADS=4

cd build

for i in {0..149}; do
    if [ ! -f "../output/expected_steps_${i}.csv.gz" ]; then
        ./yam "$i" 7&
    fi
done

wait
cd ..
./combine.sh
