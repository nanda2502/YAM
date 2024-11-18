cd output
FILENAME="expected_steps_7.csv.gz"
STRATEGY="Conformity"
zcat "$FILENAME" | grep -v "$STRATEGY" | gzip > "$FILENAME"