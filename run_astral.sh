#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=astral4
#SBATCH --mem=100G
#SBATCH --time=150:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

cd $SLURM_SUBMIT_DIR

module load ASTRAL/5.6.1-Java-1.8.0_144

#java -jar /apps/eb/ASTRAL/5.6.1-Java-1.8.0_144/astral.5.6.1.jar --input GeneTrees.tre --output Astral.out.tre

#run bootstrapping on astral tree
#java -jar /apps/eb/ASTRAL/5.6.1-Java-1.8.0_144/astral.5.6.1.jar --input Autosome_GeneTrees.tre -b Autosome_RAxML_bootstrap_filepath.txt --output Autosome_Astral.bootstrap.out.tre
#java -jar /apps/eb/ASTRAL/5.6.1-Java-1.8.0_144/astral.5.6.1.jar --input Astral.out.tre -b bins_inclusion_bootstraps.txt --output Astral_bins_inclusion_bootstrap.out

ml Python/2.7.18-GCCcore-10.2.0
#python astral_sensitivity.py
python astral_sensitivity_edited.py
