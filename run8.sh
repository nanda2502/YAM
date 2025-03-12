#!/bin/bash
#SBATCH -p genoa
#SBATCH --array=1-7
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=7     # 7 tasks per node
#SBATCH --cpus-per-task=27      # 27 CPUs per task (7×27=189, close to 192 cores)
#SBATCH --mem-per-cpu=1792M     # ~48G per task (336G ÷ 7 tasks = 48G per task, 48G ÷ 27 CPUs ≈ 1.8G per CPU)
#SBATCH -t 02:00:00
#SBATCH --output=slurm-%A.out   # Avoid having multiple output files

export OMP_NUM_THREADS=27  # Match cpus-per-task

cd build

total_tasks=49  # 0 to 48
tasks_per_group=7  # 49/7 rounded
start=$(((SLURM_ARRAY_TASK_ID - 1) * tasks_per_group ))
end=$(( SLURM_ARRAY_TASK_ID * tasks_per_group - 1))

# Ensure we don't exceed total tasks
if [ $end -gt 48 ]; then
    end=48
fi

# Create a file to track running processes
running_pids="/tmp/running_pids_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
touch $running_pids

# Launch all tasks for this array job
for i in $(seq $start $end); do
    if [ ! -f "../output/expected_steps_${i}.csv.gz" ]; then
        ./yam "$i" 50 &
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