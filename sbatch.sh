#!/bin/bash -l
#SBATCH  -N1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH  -J mea_ss
#SBATCH   -q regular
#SBATCH  --gpus-per-task=1
#SBATCH  --cpus-per-task=64 -C gpu -A  m2043_g
#SBATCH --image=rohanmalige/rohan_si-98:v2

# OpenMP settings:
export OMP_NUM_THREADS=24
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export KILOSORT2_PATH="/pscratch/sd/m/mpatil1/Kilosort2"
export HDF5_PLUGIN_PATH='/pscratch/sd/m/mpatil1/hdf5_plugin'
#kilo sort export and hdf5 export paths

# Run the application:
time srun  --output=$SCRATCH/output_log_files/%A.out --error=$SCRATCH/output_log_files/%A.err shifter python3 /pscratch/sd/m/mpatil1/MEA_Analysis_V1/IPNAnalysis/mea_analysis_pipeline.py /pscratch/sd/m/mpatil1/Data/CHD8_2/CHD8_2/230915/  -r /pscratch/sd/m/mpatil1/CHD82_ref_seventhtrial.xlsx