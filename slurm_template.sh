#!/bin/bash
#SBATCH --job-name=1j8h_1
#SBATCH -p multi
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --gres=gpu:1

# Decide the software version
export LD_LIBRARY_PATH=/public/software/lib/:$LD_LIBRARY_PATH
source /public/software/profile.d/apps_gromacs_2023.2.sh

echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"