#!/bin/bash
#SBATCH -A m2043                      # Account
#SBATCH -q regular                      # Queue
#SBATCH --constraint=cpu               # CPU only job
#SBATCH -t 7:30:00                      # Time limit
#SBATCH -N 1                            # Number of nodes
#SBATCH --mail-type=ALL                 # Send email on all job events
#SBATCH -o /pscratch/sd/a/adammwea/workspace/RBS_axonal_reconstructions/pipeline_scripts/logs/%x-%j.out
#SBATCH -e /pscratch/sd/a/adammwea/workspace/RBS_axonal_reconstructions/pipeline_scripts/logs/%x-%j.err
#SBATCH --image=adammwea/axonkilo_docker:v7

# Create logs directory if it doesn't exist
mkdir -p "$(dirname "$0")/logs"

if [ -z "$SLURM_JOB_ID" ]; then
    # Not running as an sbatch job
    shifter --image=adammwea/axonkilo_docker:v7 python /pscratch/sd/a/adammwea/workspace/RBS_axonal_reconstructions/pipeline_scripts/run_pipeline_login_latest.py > "$(dirname "$0")/logs/run_pipeline_login_latest.log" 2>&1
else
    # Running as an sbatch job
    shifter --image=adammwea/axonkilo_docker:v7 python /pscratch/sd/a/adammwea/workspace/RBS_axonal_reconstructions/pipeline_scripts/run_pipeline_login_latest.py
fi