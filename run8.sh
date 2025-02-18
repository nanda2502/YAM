#!/bin/bash
#SBATCH -p genoa
#SBATCH --array=1-8
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=6
#SBATCH -t 00:45:00
#SBATCH --output=slurm-%A.out

export OMP_NUM_THREADS=6

cd build

total_tasks=2045  # 0 to 2044
tasks_per_group=256  # Roughly 2045/8 rounded up
max_concurrent_tasks=256  # 32 tasks/node * 8 nodes
start=$(((SLURM_ARRAY_TASK_ID - 1) * tasks_per_group ))
end=$(( SLURM_ARRAY_TASK_ID * tasks_per_group - 1))

# Ensure we don't exceed total tasks
if [ $end -gt 2044 ]; then
    end=2044
fi

# Create a file to track running processes
running_pids="/tmp/running_pids_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
touch $running_pids

# Launch all tasks for this array job at once since we have enough CPU capacity
for i in $(seq $start $end); do
    if [ ! -f "../output/expected_steps_${i}.csv.gz" ]; then
        ./yam "$i" 8 &
        echo $! >> $running_pids
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
