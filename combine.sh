num_nodes=$1

cd output
zcat expected_steps_0.csv.gz | head -n 1 > header.csv

# Extract data from all files, skipping headers, and append it below the extracted header
for file in expected_steps_*.csv.gz; do
    zcat "$file" | tail -n +2 >> combined_data.csv
done

cat header.csv combined_data.csv > expected_steps.csv

rm header.csv combined_data.csv expected_steps_*

mv expected_steps.csv expected_steps_${num_nodes}.csv