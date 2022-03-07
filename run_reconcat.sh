#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=concat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu
#SBATCH --time=100:00:00

cd $SLURM_SUBMIT_DIR

module load Python

python reConCat_2019_05_08.py mafft_output/
