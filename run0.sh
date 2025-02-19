#!/bin/bash
#SBATCH -p genoa
#SBATCH --nodes=1
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=6
#SBATCH -t 01:00:00          
#SBATCH --output=slurm-%A.out

export OMP_NUM_THREADS=6

cd build
if [ ! -f "../output/expected_steps_0.csv.gz" ]; then
    ./yam "0" $1
fi