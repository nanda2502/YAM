cd output

csv_count=$(ls -1 yam_out_*.csv 2>/dev/null | wc -l)

suffix="${1:-}"

if [ -z "$suffix" ]; then
   for n in {3..100}; do
       if [ -f "../data/adj_mat_${n}.csv" ]; then
           line_count=$(wc -l < "../data/adj_mat_${n}.csv")
           if [ "$csv_count" -eq "$line_count" ]; then
               suffix=$n
               break
           fi
       fi
   done
fi

head -n 1 yam_out_0.csv > header.csv

for file in yam_out_*.csv; do
   tail -n +2 "$file" >> combined_data.csv
done

cat header.csv combined_data.csv > yam_out.csv

rm header.csv combined_data.csv

if [ -n "$suffix" ]; then
   mv yam_out.csv yam_out_${suffix}.csv
fi