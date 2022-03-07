#this script will rank the loci based on euclidean distance and create a file that has their name in order

import os

loci_dict = {}
out_file = open("bins2345678910_genenames.txt", 'w')
for files in os.listdir("/scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/Try/bins2345678910_RAxML/"):
    if files.startswith("RAxML_bestTree."):
#        count+=1
        name_split = files.split(".")
        loci_name = name_split[1].split("_")[0]
#        loci_dict[count] = loci_name
        print >>out_file, loci_name
out_file.close()

#dist_list = []
#for files in os.listdir("/scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/Try/"):
#    if files.endswith(".csv"):
#        eucl_file = open("/scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/Try/"+files, 'r')
#        for line in eucl_file:
#            things = line.split(',')
#            number = things[0][1:2]
#            eucl_dist = things[1][1:-1]
#            print number, eucl_dist

