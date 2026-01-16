#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks=4               # 4 Tasks = 1 Task per GPU
#SBATCH --time=01:30:00          # 1.5 Hours is usually enough for 4-well parallel
#SBATCH -J mea_fast_run
#SBATCH -q regular
#SBATCH --gpus-per-node=4        # Request full node (all 4 GPUs)
#SBATCH -C gpu
#SBATCH -A m2043_g
#SBATCH --image=mandarmp/benshalomlab_spikesorter:latest
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# --- 1. Load Parallel Utility ---
#module load parallel

# --- 2. Environment Setup ---
export OMP_NUM_THREADS=8         # 32 cores / 4 tasks = 8 threads per task
export HDF5_PLUGIN_PATH='/pscratch/sd/m/mpatil1/hdf5_plugin'
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"

# --- 3. Inputs (Pass these as arguments) ---
# Usage: sbatch submit_per_run.slurm /path/to/recording_folder
RUN_DIR="$1"
OUTPUT_ROOT="/pscratch/sd/m/mpatil1/MEA_Analysis/AnalyzedData"
WORKER_SCRIPT="/pscratch/sd/m/mpatil1/MEA_Analysis/IPNAnalysis/mea_analysis_routine.py"

# Auto-detect the .raw.h5 file
H5_FILE=$(find "$RUN_DIR" -name "*.raw.h5" | head -n 1)

if [ -z "$H5_FILE" ]; then
    echo "CRITICAL: No .raw.h5 file found in $RUN_DIR"
    exit 1
fi

echo "Target Run: $H5_FILE"

# --- 4. Define Worker Function ---
# This function runs a SINGLE well. GNU Parallel will spawn 4 of these at once.
process_well() {
    local well_id=$1
    echo "[$(date)] Starting $well_id on a GPU..."
    
    # --gpus-per-task=1 ensures each Python process gets its OWN unique GPU
    srun --exact --gpus-per-task=1 -n 1 --cpus-per-task=8 shifter python3 "$WORKER_SCRIPT" \
        "$H5_FILE" \
        --well "$well_id" \
        --output-dir "$OUTPUT_ROOT" \
        --checkpoint-dir "$OUTPUT_ROOT/checkpoints" \
        --sorter "kilosort4" \
        --clean-up \
        --export-to-phy \
        --force-restart

    echo "[$(date)] Finished $well_id"
}

export -f process_well
export WORKER_SCRIPT H5_FILE OUTPUT_ROOT

# --- 5. Determine Wells to Process ---
# Standard MaxTwo plates usually use well000 to well005. 
# We list them all; the pipeline will simply fail fast/skip if a well is empty.
WELLS=("well000" "well001" "well002" "well003" "well004" "well005")

# --- 6. Execute Parallel ---
# -j 4 : Run 4 jobs simultaneously (matching the 4 GPUs)
# --delay 2 : Wait 2s between starts to prevent file read collisions
echo "Launching parallel analysis on 4 GPUs..."
parallel -j 4 --delay 2 process_well ::: "${WELLS[@]}"

echo "Run Complete."