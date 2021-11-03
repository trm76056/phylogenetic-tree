#!/bin/bash                                                                                                                           
#SBATCH --job-name=sample                                                                                                            
#SBATCH --partition=highmem_30d_p                                                                                                     
#SBATCH --mail-type=ALL                                                                                                               
#SBATCH --mail-user=trm76056@uga.edu                                                                                                  
#SBATCH --ntasks=1                                                                                                                    
#SBATCH --cpus-per-task=8                                                                                                             
#SBATCH --mem=500G                                                                                                                    
#SBATCH --time=500:00:00                                                                                                              
#SBATCH --output=orp.%j.out                                                                                                           
#SBATCH --error=orp.%j.err                                                                                                            

sp="sample"

cd /scratch/trm76056/PhyloTreeData/raw_data/$sp/raw_reads

read1="sample_R1.fastq.gz"
read2="sample_R2.fastq.gz"                                                                                  

singularity run /apps/singularity-images/orp_2.2.6.sif bash -c "source activate orp &&                                                
/home/orp/Oyster_River_Protocol/oyster.mk \                                                                                           
MEM=500 \                                                                                                                             
CPU=8 \                                                                                                                               
READ1=$read1 \                                                                                                                        
READ2=$read2 \                                                                                                                        
RUNOUT=$sp "
