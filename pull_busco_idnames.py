#this script will pull out the names of the transcripts that matched to busco sequences
#will be used to pull out transcripts using samtools faidx in another script to have that as input for transdecoder

import os

for dirs in os.listdir("busco_results/"):
     things = dirs.split("_")
     print(things)
     if len(things)==3:
          sp = things[2]
     if len(things)==4:
          sp = things[2]+"_"+things[3]
     table_file = open("busco_results/"+dirs+"/full_table_busco_"+sp+".tsv", 'r')
     out_file = open(sp+"_buscoID.txt", 'w')
     for line in table_file:
          if line.startswith("#"):
               pass
          else:
               stuff = line.split()
               if len(stuff) > 2:
                    busco_id = stuff[0]
                    seq_id = stuff[2]
                    print >>out_file, seq_id