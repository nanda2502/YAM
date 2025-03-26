#!/bin/bash
#SBATCH -p genoa
#SBATCH --array=0-16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH -t 04:45:00
#SBATCH --output=slurm-%A.out

export OMP_NUM_THREADS=192

cd build

task_id=${SLURM_ARRAY_TASK_ID}

# Check if the output already exists
if [ ! -f "../output/expected_steps_${task_id}.csv.gz" ]; then
    # Run the task
    ./yam "$task_id" 2
fi