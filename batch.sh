#! /bin/bash
#SBATCH -J n-step
#SBATCH -t 00:10:00
#SBATCH -p thin
#SBATCH -n 512

EXECUTABLE="$HOME/YAM/build/yam.exe"
cd build

for n in $(seq 3 8); do
    "$EXECUTABLE" "$n" &
done