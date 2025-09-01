#!/bin/bash

min_p=2
max_p=8

min_t=1
max_t=8

declare -a sizes=("8GB")
filesize=8GB
declare -a trim_modes=("3" "5")
trim_mode=5
declare -a sequences=("TATA" "AAATAAATATAGAACC")
sequence=TATA

echo "Process counts vertically from $min_p to $max_p"
echo "Thread counts horizontally from $min_t to $max_t"

for filesize in "${sizes[@]}"; do
echo ""
echo "===== NEW SIZE: $filesize  ====="
# TODO: Compress FASTQ if NGSC doesn't exist

for trim_mode in "${trim_modes[@]}"; do
echo ""
echo "--- Trimming $trim_mode' end ---"
echo ""

for sequence in "${sequences[@]}"; do
echo "Adapt-seq: $sequence"

    for (( p=min_p; p <= max_p; p*=2 ))
    do
        echo -n "($p)    "

        for (( t=min_t; t <= $max_t; t*=2 ))
        do
            # echo -n "mpiexec -np $p ./main -i $t out${filesize}.ngsc testing/trimmed${trim_mode}_${filesize}.fastq -trim -$trim_mode $sequence -debug"
            # echo ""
            output=$(mpiexec -np $p ./main -i $t out${filesize}.ngsc testing/trimmed${trim_mode}_${filesize}.fastq -trim -$trim_mode $sequence -debug)
            output=$((echo "$output" | grep -P "Max runtime: ") | grep -oP "\d+(\.\d+)?")
            echo -n "$output, "
        done

        echo ""
    done

    echo ""
done
done
done