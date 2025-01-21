export OMP_NUM_THREADS=40

num_nodes=$1

file="./data/adj_mat_${num_nodes}.csv"

max_i=$(($(wc -l < "$file") - 1))

cd build

for i in $(seq 0 $max_i); do
    if [ ! -f "../output/expected_steps_${i}.csv.gz" ]; then
        ./yam "$i" $num_nodes&
    fi
done

wait 
cd .. 
./combine.sh $num_nodes