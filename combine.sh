cd output

gz_count=$(ls -1 expected_steps_*.csv.gz 2>/dev/null | wc -l)

suffix="${1:-}"

if [ -z "$suffix" ]; then
   for n in {3..100}; do
       if [ -f "../data/adj_mat_${n}.csv" ]; then
           line_count=$(wc -l < "../data/adj_mat_${n}.csv")
           if [ "$gz_count" -eq "$line_count" ]; then
               suffix=$n
               break
           fi
       fi
   done
fi

zcat expected_steps_0.csv.gz | head -n 1 > header.csv

for file in expected_steps_*.csv.gz; do
   zcat "$file" | tail -n +2 >> combined_data.csv
done

cat header.csv combined_data.csv > expected_steps.csv

rm header.csv combined_data.csv expected_steps_*.gz

if [ -n "$suffix" ]; then
   mv expected_steps.csv expected_steps_${suffix}.csv
fi