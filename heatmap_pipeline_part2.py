# Part two of the pipeline
# Written by Shannon Keating
# On python version 3.71.
# April 2019

import os 
import sys
import subprocess
import re
from natsort import natsorted

# To run this script:
# python part2.py [directory_with_tiger_output_nexus_files]
# for example
# python part2.py /shannon/data/tiger/out

print("Splitting and concatenating nexus alignments based on Tiger-provided bins")
print("\n")
print("Would you like to clean up folder once done, leaving only the final concatenated alignments? Please enter yes or no.")
cleanup = input()
print(cleanup)
if cleanup.lower() == "yes":
	print("This program will erase all intermediate files once finished.")
elif cleanup.lower() == "no":
	print("This program will keep all intermediate files once finished.")
else:
	print("Answer not understood. Quiting program")
	exit()
	
	
# get a list of all the alignment files in a directory. Check that they are nexus files
fasta_directory = sys.argv[1]
os.chdir(fasta_directory)
alignments = []
folder_contents = os.listdir(".") 
for file in folder_contents:
	if file.endswith(".nex") or file.endswith(".nexus"):
		alignments.append(file)

# align will be the names pulled from the folder, i.e. the contents of alignments
# file will be the file pulled and read in each time
# contents will be the file contents, i.e. the actual strings
for align in alignments:
	# make the partition file for each alignment
	print("Making partition file for all nexus files in folder")
	print(align)
	file = open(align, 'r')
	name_split = align.split(".")
	name = name_split[0]
	partition = open(name + "_partition.txt", 'w')
	for lines in file:
		bins = lines.strip()
		if bins.startswith('Charset'):
			partition.write(bins[8:-1]+ "\n")
	partition.close()
	
	# fix the 'Matrix' problem AMAS seems to have (won't run if the 'm' in 'matrix' is capitalized. Not sure why this is
	file = open(align, 'r')
	filedata = file.read() 
	file.close()
	newdata = filedata.replace("Matrix","matrix") # use this for AMAS
	file = open(align,'w') # open the file again to overwrite it
	file.write(newdata)
	file.close()
	
	# running AMAS split for all files
	print("Splitting alignment files based on partitions......")
	AMAS_part = "/Applications/AMAS-master/amas/AMAS.py split -f nexus-int -d dna -u nexus"
	input = align
	partition = (name + "_partition.txt")
	AMAS_call = AMAS_part+" -i {0} -l {1}".format(input,partition)
	print("AMAS call is " + AMAS_call)
	run_call = subprocess.call(AMAS_call, shell = True)


	# concatenating files back together based on bin
	print("Concatenating....")
	folder_contents = os.listdir() 
	AMAS_concat = "/Applications/AMAS-master/amas/AMAS.py concat -f nexus -d dna -u phylip"
	print("Running AMAS for " + name + "\n")
	bin_list = []
	for input in folder_contents:
		if input.startswith(name+"_Bin"):
			bin_list.append(input)
	bin_list = natsorted(bin_list)
	print(bin_list)
	to_concat = ''
	for i in bin_list:
		filename = i.split("-")
		filename = filename[0]
		to_concat = str(to_concat + " " + i).rstrip()
		#Build up the AMAS string
		AMAS_call = AMAS_concat+" -i {0} -p {1} -t {2}".format(to_concat,filename+"_partition.txt",filename+"_concatenated.phy")
		print("AMAS call is " + AMAS_call + "\n")
		print("Running AMAS")
		run_call = subprocess.call(AMAS_call, shell = True)


	# Cleaning up the folder
	if cleanup.lower() == "yes":
		print("\n\nCleaning up folder and removing individual bin files")
		files = os.listdir()
		for i in files:
			if re.search("Bin.*-out.nex", i):
				os.remove(i)

print("\nFinished pipeline part 2")