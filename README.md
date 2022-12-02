# phylogenetic-tree

This is the  overview of scripts and the pipeline used for this project.

TRANSCRIPTOME ASSEMBLY AND DATASET GENERATION

- Once reads were received, paired end RNA sequencing reads were used as input for the OysterRiverProtocol.sh script to generate a transcriptome assembly.
- We then ran the transcriptomes through the metazoa_odb9 database from BUSCO to pull out single-copy orthologs.
  - This was the run_busco.sh script
- We pulled the genome assemblies for five additional species.
  - https://web.stanford.edu/~bkim331/files/genomes/ was used for D. grimshawi, D. pruinosa, D. virilis, and D. quadrilineata
  - D. innubila genome was downloaded from NCBI (GenBank accession GCA_004354385.2)
  - These genomes were run through the metazoa_odb9 database from BUSCO to pull out the orthologs.
  - We had to create a script that would pull out the exact BUSCO sequences since it seemed to be pulling out full contigs and was impacting the TOAST pipeline.
    - pull_busco_genomes.py was used to get the sequence information
    - pull_nuc_seq.sh helped grab the full sequence from the genome for those orthologs
    - matchbusco_mafft.py appended those ortholog sequences to the extracted files from TOAST
- TOAST (cite) is an ortholog alignment program that uses the output from BUSCO to generate an alignment.
  - We used ParseBuscoResults and ExtractBuscoSeqs on the transcriptome data.
  - The extracted fasta files that included both transcriptome and genome sequence data were run through MafftOrientAlign, MissingDataTable, and SuperAlign to generate a final concatenated alignment file for all metazoan BUSCO orthologs.
- We split the dataset into X-linked and autosome orthologs.
  - We used the D. melanogaster X chromosome sequence to generate a Blast database and pulled out the ortholog sequences that blasted to the D. mel X. 
  - Orthologs that didn't match to the X chromosome were presumed to be autosome orthologs.
  - The aligned mafft files for the two ortholog datasets were used as input for TOAST's SuperAlign function to create a concatendated alignment for the two datasets.

PHYLOGENETICS

- The concatenated datasets that were generated were used as input for a maximum likelihood analysis using IQ-Tree v2.0.6.
  - We used Model Finder to determine the best fit model and performed 1000 bootstraps on the bests fit model.
  - Trees were viewed using FigTree v1.4.4.
- Individual gene trees were generated with RAxML v8.1.12 with the GTRCAT model and 100 bootstraps.
  - These individual gene trees were concatenated and used to estimate the species tree using the coalescent-based method in ASTRAL. 
- ML and coalescent trees were generated for all three datasets (whole, X, autosome). 
- We inferred anscestral states using RASP.
  - We did this separately for condensed phylogenies that consisted of either the flies that were tested for amanita tolerance or the Natural toxin mix tolerance. These were edited in Inkscape.
- We performed a sensitivity analysis to identify orthologs that best supported our species tree.
  - This was done first by using tiger (tiger.sh) to organize sites based on site disagreement.
  - Then, we used AMAS to concatenate binned sequences and sequentially added bins to add sites.
  - We regenerated ML trees for each alignment with RAxML. 
    - We then estimated pairwise distance among trees with treeCMP. Plotted these with cmdscale in R and calculated euclidean distance.
    - generated new ML trees with bootstrap support values for new alignments based on euclidean distance and binning.
  - Used bootstrapped gene trees to get support for coalescent tree in RAxML. 
  - Generated heatmaps with custom script from Dowdy et al. to get heatmap for each node of the tree. Manually edited in Inkscape.
