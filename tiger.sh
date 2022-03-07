#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=folder_9
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu
#SBATCH --mem=100G

cd $SLURM_SUBMIT_DIR

ml Python/3.8.6-GCCcore-10.2.0
#ml Python/2.7.18-GCCcore-10.2.0

#index aligned sequence using tiger
#for files in mafft_aligned/*.fasta
#do
    #index aligned sequence using tiger
    echo $files
    #python tiger.py index -i $files -o $files.tiger
    #get rate using tiger
    #python tiger.py rate -i $files.tiger.ref.ti -o $files.tiger.rate
    #run output for tiger
#    biotiger-master/tiger output -i $files.tiger.rate.gr --mask -exc 1 -fa $files -o $files.Bins2345678910.tiger.output
#done

#python reConCat_2019_05_08.py folder_
biotiger-master/tiger index -i EOG091G0EKA.fasta -o EOG091G0EKA.fasta.tiger
biotiger-master/tiger rate -i EOG091G0EKA.fasta.ref.ti -o EOG091G0EKA.fasta.tiger.rate
biotiger-master/tiger output -i EOG091G0EKA.fasta.rate.gr --mask -exc 1,3,4,5,6,7,8,9,10 -fa EOG091G0EKA.fasta -o EOG091G0EKA.fasta.Bins2.tiger.output
