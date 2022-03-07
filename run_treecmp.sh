#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=treecmp2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=100:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu

cd $SLURM_SUBMIT_DIR

ml Java/1.8.0_144

java -jar ~/TreeCmp/bin/TreeCmp.jar -m -d tt -P -i bins2345678910_RAxML.nexus -o bins2345678910_treecmp.out
