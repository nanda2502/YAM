#!/bin/bash
#SBATCH -p genoa
#SBATCH --array=0-255           # Process 256 array elements
#SBATCH --nodes=1               # Each array job uses 1 node
#SBATCH --ntasks-per-node=8     # 8 tasks per node
#SBATCH --cpus-per-task=24      # 24 CPUs per task
#SBATCH -t 05:00:00
#SBATCH --output=slurm-%A.out   # Single output file for all array tasks

export OMP_NUM_THREADS=24  # Match cpus-per-task

cd build

total_tasks=2045            
tasks_per_array=8

# Calculate start and end indices for this array task
start=$((SLURM_ARRAY_TASK_ID * tasks_per_array))
end=$((start + tasks_per_array - 1))

# Ensure we don't exceed total tasks
if [ $end -ge $total_tasks ]; then
    end=$((total_tasks - 1))
fi

echo "Array task ${SLURM_ARRAY_TASK_ID} processing indices ${start} to ${end}"

# Create a file to track running processes
running_pids="/tmp/running_pids_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
touch $running_pids

# Launch all tasks for this array job
for i in $(seq $start $end); do
    if [ ! -f "../output/expected_steps_${i}.csv.gz" ]; then
        echo "Starting task $i"
        ./yam "$i" 8 &
        echo $! >> $running_pids
    else
        echo "Skipping task $i (output already exists)"
    fi
done

# Wait for all processes to complete
while read -r pid; do
    if kill -0 $pid 2>/dev/null; then
        wait $pid
    fi
done < $running_pids

# Clean up
rm $running_pids