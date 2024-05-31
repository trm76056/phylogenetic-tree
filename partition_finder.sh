#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=partitionfinder
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=100:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu

cd $SLURM_SUBMIT_DIR

ml PartitionFinder/2.1.1-foss-2019b-Python-2.7.16

time python $EBROOTPARTITIONFINDER/PartitionFinder.py -p 4 --RAxML --min-subset-size 5000