num_nodes=$1

file="./data/adj_mat_${num_nodes}.csv"

max_i=$(($(wc -l < "$file") - 1))

missing=0

for i in $(seq 0 $max_i); do
    if [ ! -f "./output/expected_steps_${i}.csv.gz" ]; then
        missing=$((missing + 1))
    fi
done

run_script="run${num_nodes}.sh"

sed -i "s/^#SBATCH -n .*/#SBATCH -n $missing/" "$run_script"

sbatch "$run_script"