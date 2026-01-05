#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --time=01:30:00          # Adjusted time (KS4 is faster than KS2)
#SBATCH -J mea_pipeline
#SBATCH -q regular               # or 'debug' for testing
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=32        # 64 is overkill for single well, 32 is plenty
#SBATCH -C gpu
#SBATCH -A m2043_g
#SBATCH --image=mandarmp/benshalomlab_spikesorter:latest
#SBATCH --output=%x_%j.out        # Saves to simpler log names (mea_pipeline_12345.out)
#SBATCH --error=%x_%j.err

# --- 1. Environment Setup ---
# Maximize CPU threading for SpikeInterface preprocessing
export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Essential for Maxwell .h5 files
export HDF5_PLUGIN_PATH='/pscratch/sd/m/mpatil1/hdf5_plugin'

# CRITICAL for Kilosort4 on PyTorch to prevent fragmentation OOM
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"

# --- 2. Input Variables (Edit These) ---
# Define the raw file and output location
INPUT_PATH="/pscratch/sd/m/mpatil1/Data/CDKL5_T1/240607/M07039/Network/000088/data.raw.h5"
OUTPUT_ROOT="/pscratch/sd/m/mpatil1/MEA_Analysis/AnalyzedData"
SCRIPT_PATH="/pscratch/sd/m/mpatil1/MEA_Analysis/IPNAnalysis/run_pipeline_driver.py"



# --- 3. Run the Pipeline ---
echo "Starting Analysis for $WELL_ID on $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"

time srun shifter python3 "$SCRIPT_PATH" \
    "$INPUT_PATH" \
    --output-dir "$OUTPUT_ROOT" \
    --checkpoint-dir "$OUTPUT_ROOT/checkpoints" \
    --sorter "kilosort4" \
    --clean-up \
    --export-to-phy

echo "Job Complete"