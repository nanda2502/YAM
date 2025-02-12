#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <num_nodes>"
    exit 1
fi

# Set expected total based on num_nodes
case "$1" in
    8) total=554 ;;
    7) total=150 ;;
    6) total=57 ;;
    *) echo "Invalid num_nodes. Use 6, 7, or 8."; exit 1 ;;
esac

# Function to draw progress bar
draw_progress_bar() {
    local width=50
    local percentage=$1
    local completed=$(( (width * percentage) / 100 ))
    local remaining=$((width - completed))
    
    printf "["
    printf "%${completed}s" | tr ' ' '#'
    printf "%${remaining}s" | tr ' ' '-'
    printf "] %3d%%" "$percentage"
}

# Function to calculate rate based on file timestamps
calculate_rate() {
    # Get list of files sorted by modification time
    files=($(ls -t ./output/expected_steps_*.csv.gz 2>/dev/null))
    count=${#files[@]}
    
    if [ $count -lt 2 ]; then
        echo "0"
        return
    fi
    
    # Get timestamps of newest and oldest files
    newest_time=$(stat -c %Y "${files[0]}")
    oldest_time=$(stat -c %Y "${files[-1]}")
    time_diff=$((newest_time - oldest_time))
    
    if [ $time_diff -eq 0 ]; then
        echo "0"
        return
    fi
    
    # Calculate files per minute
    rate=$(echo "scale=2; $count * 60 / $time_diff" | bc)
    echo "$rate"
}

while true; do
    # Count current files
    current_count=$(ls -1 ./output/expected_steps_*.csv.gz 2>/dev/null | wc -l)
    
    # Calculate percentage
    if [ $total -gt 0 ]; then
        percentage=$((current_count * 100 / total))
    else
        percentage=0
    fi
    
    if [ $current_count -gt 0 ]; then
        # Calculate rate based on file timestamps
        rate=$(calculate_rate)
        
        if [ "$rate" != "0" ]; then
            # Calculate remaining time
            remaining_files=$((total - current_count))
            mins_remaining=$(echo "scale=0; $remaining_files / $rate" | bc)
            
            printf "\033[K"
            draw_progress_bar $percentage
            printf " | %d/%d files | %.1f files/min | ETA %d mins\r" \
                   "$current_count" "$total" "$rate" "$mins_remaining"
        else
            printf "\033[K"
            draw_progress_bar $percentage
            printf " | %d/%d files | Calculating rate...\r" \
                   "$current_count" "$total"
        fi
    else
        printf "\033[K"
        draw_progress_bar $percentage
        printf " | Waiting for first file...\r"
    fi
    
    # Exit if we're done
    if [ $current_count -eq $total ]; then
        echo -e "\nComplete! All files processed."
        newest_file=$(ls -t ./output/expected_steps_*.csv.gz 2>/dev/null | head -1)
        oldest_file=$(ls -t ./output/expected_steps_*.csv.gz 2>/dev/null | tail -1)
        start_time=$(stat -c %Y "$oldest_file")
        end_time=$(stat -c %Y "$newest_file")
        total_time=$((end_time - start_time))
        echo "Total time: $((total_time / 60)) minutes $((total_time % 60)) seconds"
        exit 0
    fi
    
    sleep 10
done