#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=mafft
#SBATCH --mem=100G
#SBATCH --ntasks=10
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu
#SBATCH --time=100:00:00

cd $SLURM_SUBMIT_DIR

ml MAFFT/7.453-GCC-8.3.0-with-extensions

mafft --adjustdirection --thread 10 --quiet ", extract_dir, "/", file_to_use, " > ", mafft_dir, "/", file_to_use