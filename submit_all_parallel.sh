#!/bin/bash

# Usage: ./submit_all.sh /path/to/root/data/folder

ROOT_DIR="$1"

if [ -z "$ROOT_DIR" ]; then
    echo "Usage: $0 <path_to_root_data_directory>"
    exit 1
fi

echo "=========================================================="
echo "SCANNING: $ROOT_DIR"
echo "Targeting only: */Network/*/data.raw.h5"
echo "=========================================================="

# 1. Find files matching the specific pattern */Network/*/data.raw.h5
find "$ROOT_DIR" -path "*/Network/*" -name "data.raw.h5" | sort | while read -r FILEPATH; do
    
    # Get the directory containing the file
    RUN_DIR=$(dirname "$FILEPATH")
    
    # Create a cleaner job name for Slurm (e.g., M05506_Network_000014)
    # Extracts the last 3 parts of the path for identification
    JOB_ID_STR=$(echo "$RUN_DIR" | awk -F/ '{print $(NF-2)"_"$(NF-1)"_"$(NF)}')
    
    echo "Found Run: $RUN_DIR"
    
    # 2. Submit the job
    # -J sets the job name so you can identify it in the queue
    sbatch -J "mea_${JOB_ID_STR}" \
           sbatch_parallel.sh "$RUN_DIR"
           
    echo "Submitted."
    echo "----------------------------------------------------------"

done

echo "Submission loop complete. Check queue with 'squeue --me'"