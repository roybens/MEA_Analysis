#!/bin/bash
#SBATCH -A m2043_g                      # Account
#SBATCH -C gpu                          # Request GPU nodes
#SBATCH -q regular                      # Queue
#SBATCH -t 6:00:00                      # Time limit
#SBATCH -N 1                            # Number of nodes
#SBATCH --gpus-per-task=1               # One GPU per task
#SBATCH --mail-type=ALL                 # Send email on all job events
#SBATCH -o /pscratch/sd/a/adammwea/workspace/RBS_axonal_reconstructions/pipeline_scripts/logs/%x-%j.out
#SBATCH -e /pscratch/sd/a/adammwea/workspace/RBS_axonal_reconstructions/pipeline_scripts/logs/%x-%j.err
#SBATCH --image=adammwea/axonkilo_docker:v7

#target well 0 for testing
shifter --image=adammwea/axonkilo_docker:v7 python /pscratch/sd/a/adammwea/workspace/RBS_axonal_reconstructions/pipeline_scripts/run_pipeline_HPC_latest.py