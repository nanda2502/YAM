export OMP_NUM_THREADS=8

cd build

for i in {0..6}; do
    if [ ! -f "../output/expected_steps_${i}.csv.gz" ]; then
        ./yam "$i" 4&
    fi
done