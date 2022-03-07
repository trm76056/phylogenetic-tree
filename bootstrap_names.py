out_file = open('Autosome_RAxML_bootstrap_filepath.txt', 'w')

import os

for files in os.listdir("/scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/Try/Autosome_RAxML_genetrees_boot/"):
     things = files.split(".")
     busco = things[1]
     for x_files in os.listdir("/scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/Try/RAxML_genetrees_boot/"):
          if x_files.startswith("RAxML_bootstrap."):
               stuff = x_files.split(".")
               if stuff[1]==busco:
                    new_file=x_files.strip()
                    print(new_file)
                    line_name = "RAxML_genetrees_boot/"+new_file
                    print >>out_file, line_name