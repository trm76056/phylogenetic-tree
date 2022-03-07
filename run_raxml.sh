#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=raxmlbipart
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu
#SBATCH --mem=10G

cd /scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/Try

#ml RAxML/8.2.12-intel-2019b-hybrid-avx2

#for x in *.phy;
#do
#     echo $x
#     raxmlHPC -s $x -n $x.out -T 4 -f a -m GTRGAMMA -x 12345 -# 100 -p 12345 -w /scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/Try/combination_sets/RAxML_boot
#     raxmlHPC -s $x -n $x.out -p 12345 -T 4 -# 100 -m GTRCAT -f a -w /scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/Try/RAxML_genetrees
#done

#raxmlHPC -s bins23456789/EOG091G0803_Bins23456789_concatenated.phy -n EOG091G0803_Bins23456789_concatenated.phy.out -T 4 -f a -m GTRGAMMA -x 12345 -# 100 -p 12345 -w /scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/Try


ml RAxML/8.2.12-foss-2019b-pthreads-avx
ml Python/2.7.18-GCCcore-10.2.0
#ml RAxML-NG/0.9.0-GCC-8.3.0
#raxmlHPC -f e -t Astral.out.tre -m GTRGAMMA -s superalign.fa -n bins_astral_bipartitions_fb.out -T 4
#raxmlHPC -f b -t RAxML_result.bins_astral_bipartitions_fb.out -z bins_astral_bs.txt -T 4 -m GTRGAMMA -n bins_astral_bipartition.out
#raxmlHPC -t Phylo.treefile -n Phylo.treefile.sens.out -p 12345 -T 4 -# 100 -m GTRCAT -f b -z bins2345678910_inclusion1_bootstrapout.txt
python run_raxml_sensitivity.py
#raxml-ng --all --msa superalign.fa --model GTR+G
