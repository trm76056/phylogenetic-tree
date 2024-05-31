#this script will be a wrapper script to create the gene tree file that is needed to run astral as well as pull the bootstrap files for those loci that are included to create the bootstrap file 

import os
import sys
import subprocess
import re
file_list_1 = []
file_list_2 = []
file_list_3 = []
file_list_4 = []
file_list_5 = []
file_list_6 = []
file_list_7 = []
file_list_8 = []
eucl_file = open("bins23_euclidean_dist.csv", 'r')
count = 0
for lines in eucl_file:
    if lines.startswith('"",'):
        pass
    else:
        things = lines.split(",")
        count +=1
        gene_name = things[1]
        gene_name = gene_name[1:-1]
#        print gene_name                                                                                                                                                                                           
        if int(count) <= 61:
            file_list_1.append(gene_name)
        if int(count) <= 122:
            file_list_2.append(gene_name)
        if int(count) <= 183:
            file_list_3.append(gene_name)
        if int(count) <= 244:
            file_list_4.append(gene_name)
        if int(count) <= 305:
            file_list_5.append(gene_name)
        if int(count) <= 366:
            file_list_6.append(gene_name)
        if int(count) <= 427:
            file_list_7.append(gene_name)
        if int(count) <= 489:
            file_list_8.append(gene_name)
eucl_file.close()

astral_path = "/apps/eb/ASTRAL/5.6.1-Java-1.8.0_144/astral.5.6.1.jar"

for x in range(1:9):
	#create bootstrap file that's needed for astral bootstrap as well as make the file needed for the astral input
	file_1 = "bins23_inclusion"+x+"_genetrees.tre"
	file_2 = "bins23_inclusion"+x+"_bootstrap.txt"
	inclusion = open(file_1,'w')
	bootstrap = open(file_2,'w')
	file_list = "file_list_"+x
	for loci in file_list:
		loci_file = open("bins23_RAxML/RAxML_bestTree."+loci+"_Bins23_concatenated.phy.out", 'r')
		for line in loci_file:
			print >>inclusion, line
		bs_file = "bins23_RAxML/RAxML_bootstrap."+loci+"Bins23_concatenated.phy.out"
		print >>bootstrap, bs_file
	inclusion.close()
	bootstrap.close()
	#run astral process
	astral_out = "bins23_inclusion"+x+"_astral.tre"
	astral_call = "java -jar "+astral_path+"--input "+file_1+" --output "+astral_out	
	print "Astral call is "+astral_call+"\n"
	run_call = subprocess.call(astral_call, shell=True)
	#run bootstrapping on astral
	boot_out = "bins23_inclusion"+x+"_astral_bs.tre"
	boot_call = "java -jar "+astral_path+"--input "+file_1" -b file_2 --output "+boot_out
	print "Bootstrapping Astral \n"
	run_call = subprocess.call(boot_call, shell=True)

