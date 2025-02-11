#!/bin/bash

# Enable job control
set -m

# Cleanup function to kill all background processes
cleanup() {
    echo "Cleaning up..."
    # Kill all background processes in the same process group
    kill $(jobs -p) 2>/dev/null
    exit 1
}

# Trap Ctrl+C (SIGINT) and call cleanup
trap cleanup SIGINT SIGTERM

export OMP_NUM_THREADS=1

num_nodes=$1
file="./data/adj_mat_${num_nodes}.csv"
max_i=$(($(wc -l < "$file") - 1))
max_parallel=$(($(nproc) - 1))

cd build || exit 1

active_jobs=0

for i in $(seq 0 $max_i); do
    if [ ! -f "../output/expected_steps_${i}.csv.gz" ]; then
        while [ $active_jobs -ge $max_parallel ]; do
            wait -n
            active_jobs=$((active_jobs - 1))
        done
        
        ./yam "$i" $num_nodes &
        active_jobs=$((active_jobs + 1))
    fi
done

wait

cd ..
./combine.sh