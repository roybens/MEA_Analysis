#!/bin/bash -l
#SBATCH  -N1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH  -J connections
#SBATCH   -q regular
#SBATCH  --cpus-per-task=128 -C cpu 
#SBATCH --image=balewski/ubu20-neuron8:v5

# OpenMP settings:
export OMP_NUM_THREADS=24
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Directories for logs
OUTPUT_DIR=$SCRATCH/MEA_Analysis/Connections/output_log_files
mkdir -p $OUTPUT_DIR
INPUT_DATA_PATH=$SCRATCH/MEA_Analysis/Connections/data/CDKL5_T1/DIV20/M08018/well5/spike_times.npy
# Run the application
time srun \
    --output=${OUTPUT_DIR}/%A.out \
    --error=${OUTPUT_DIR}/%A.err \
    shifter python3 /pscratch/sd/m/mpatil1/MEA_Analysis/Connections/UoI_Lasso4.py \
    --data_path ${INPUT_DATA_PATH}