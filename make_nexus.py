import os

output = open("bins2345678910_RAxML.nexus", 'w')

count = 0
for files in os.listdir("bins2345678910_RAxML/"):
    if files.startswith("RAxML_bestTree."):
        count +=1
        name = files.split(".")
        gene_name = name[1].split("_")[0]
        treefile = open("bins2345678910_RAxML/"+files, 'r')
        for lines in treefile:
            newlines = "tree_"+gene_name+" = "+lines.strip()
            print>>output, newlines
        treefile.close()
output.close()
