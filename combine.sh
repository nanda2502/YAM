cd output

gz_count=$(ls -1 expected_steps_*.csv.gz 2>/dev/null | wc -l)

# Check each possible adj_mat file and compare line counts
for n in {3..8}; do
    if [ -f "../data/adj_mat_${n}.csv" ]; then
        line_count=$(wc -l < "../data/adj_mat_${n}.csv")
        if [ "$gz_count" -eq "$line_count" ]; then
            num_nodes=$n
            break
        fi
    fi
done

if [ -z "$num_nodes" ]; then
    echo "Warning: Could not determine num_nodes - using default output name"
    # Set num_nodes to empty string so we use the default name
    num_nodes=""
fi

zcat expected_steps_0.csv.gz | head -n 1 > header.csv

# Extract data from all files, skipping headers, and append it below the extracted header
for file in expected_steps_*.csv.gz; do
    zcat "$file" | tail -n +2 >> combined_data.csv
done

cat header.csv combined_data.csv > expected_steps.csv

rm header.csv combined_data.csv expected_steps_*

# Only append num_nodes to the filename if it's available
if [ -n "$num_nodes" ]; then
    mv expected_steps.csv expected_steps_${num_nodes}.csv
fi