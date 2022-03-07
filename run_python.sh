#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=python
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=5:00:00

cd $SLURM_SUBMIT_DIR

ml Python/2.7.18-GCCcore-10.2.0

python get_individual_genetrees.py
