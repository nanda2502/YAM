#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <num_nodes>"
    exit 1
fi


num_nodes="$1"


max_index=$(( $(wc -l < "./data/adj_mat_${num_nodes}.csv") - 1 ))

# Initialize counter for missing files
missing_count=0

# Loop through expected files and count missing ones
for i in $(seq 0 $max_index); do
    if [ ! -f "./output/expected_steps_${i}.csv.gz" ]; then
        echo "Missing file: expected_steps_${i}.csv.gz"
        missing_count=$((missing_count + 1))
    fi
done

# Report the number of missing files
echo $missing_count