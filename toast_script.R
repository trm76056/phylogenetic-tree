rm(list=ls())

td <- "~/Documents/Research/My_TOAST_Project/new_toast"
fd <- "~/Documents/Research/My_TOAST_Project/new_toast/fasta"
bs <- "~/Documents/Research/My_TOAST_Project/run_BUSCO.py"
bd <- "~/Documents/Research/My_TOAST_Project/new_toast/busco_results"
ed <- "~/DocumentsResearch/My_TOAST_Project/new_toast/extracted"
md <- "~/Documents/Research/My_TOAST_Project/new_toast/mafft_aligned"
od <- "~/Documents/Research/My_TOAST_Project/new_toast/laurasiatheria_odb9"
ad <- "~/Documents/Research/My_TOAST_Project/new_toast/mafft_aligned"
go_oro <- "~/Documents/Research/My_TOAST_Project/goslim_agr.obo"
ortho <- "~/Documents/Research/My_TOAST_Project/new_toast/laurasiatheria_odb9/info/laurasiatheria_314145_OrthoDB9_orthogroup_info.txt"
cpu <- 12

setwd(td)

#### Download the sequences ####

EntrezDownload <- function(txid , fasta_dir, minimumSeq = 350, maximumSeq = NULL){
  
  pull1_again <- function(){
    Sys.sleep(0.5)
    tryCatch(pull1(), error = function(e) {pull1_again()})
  }
  
  pull1 <- function(){
    write_me <- entrez_fetch(db = "nucleotide", web_history = cookie, rettype = "fasta", retstart = j, query_key = qk, retmax = 250)
  }
  
  pull2_again <- function(){
    Sys.sleep(0.5)
    tryCatch(pull2(), error = function(e) {pull2_again()})
  }
  
  pull2 <- function(){
    write_me <- entrez_fetch(db = "nucleotide", web_history = cookie, rettype = "fasta", retstart = j, query_key = qk, retmax = dif)
  }
  
  cat("STEP 1 - GATHERING DATA FROM NCBI\n")
  txid_int <- paste0("txid", txid, "[Organism:exp]") #txid comes from sourcing config.R
  taxon_search <- entrez_search(db = "taxonomy", term = txid_int, retmax = 31000) #look into retmax
  cat("The search term", txid_int, "returned", taxon_search$count, "species hits
        DOWNLOADING SEQUENCES FOR ANY SPECIES WITH >", minimumSeq, "SEQUENCES\n")
  
  for (i in 1:length(taxon_search$ids)) {
    building_search_term <- paste("txid", taxon_search$ids[i], sep = "")
    search_term <- paste(building_search_term, " AND biomol mrna[prop]", sep = "")
    nuc_search <- entrez_search(db = "nucleotide", term = search_term, use_history = TRUE) #look into the retmax
    cookie <- nuc_search$web_history
    qk <- cookie$QueryKey
    
    if (is.null(maximumSeq) == FALSE && (nuc_search$count > minimumSeq)) {
      dif <- maximumSeq %% 250
      taxize_summ <- entrez_summary(db = "taxonomy", id = taxon_search$ids[i])
      genus <- taxize_summ$genus
      species <- taxize_summ$species
      file_name <- paste(genus, species, sep = "_")
      fasta_file_name <- paste(file_name, "fasta", sep = ".")
      cat("\nDownloading", maximumSeq,"sequences for", file_name,
          "species number", i, "of", taxon_search$count, "\n")
      pb1b <- txtProgressBar(min = 0, max = maximumSeq, style = 3)
      
      for (j in seq(from=1, to=maximumSeq, by=250)) { #retmax 500 downloads faster but times out occaisionally
        j <<- j
        tryCatch(write_me <- entrez_fetch(db = "nucleotide", web_history = cookie, rettype = "fasta", retstart = j, query_key = qk, retmax = 250), error = function(e) {pull1_again()})
        write(write_me, file = paste0(fasta_dir, "/", fasta_file_name), append=TRUE)
        setTxtProgressBar(pb1b, j)
      }
      close(pb1b)
      
    } else {
      
      if (nuc_search$count > minimumSeq){ #probably change this number to make sure only things with a large number of sequences are included
        dif <- nuc_search$count %% 250
        taxize_summ <- entrez_summary(db = "taxonomy", id = taxon_search$ids[i])
        genus <- taxize_summ$genus
        species <- taxize_summ$species
        file_name <- paste(genus, species, sep = "_")
        fasta_file_name <- paste(file_name, "fasta", sep = ".")
        cat("\nDownloading", nuc_search$count,"sequences for", file_name,
            "species number", i, "of", taxon_search$count, "\n")
        pb1b <- txtProgressBar(min = 0, max = nuc_search$count, style = 3)
        
        for (j in seq(from=1, to=nuc_search$count-250, by=250)) { #retmax 500 downloads faster but times out occaisionally
          j <<- j
          tryCatch(write_me <- entrez_fetch(db = "nucleotide", web_history = cookie, rettype = "fasta", retstart = j, query_key = qk, retmax = 250), error = function(e) {pull1_again()})
          write(write_me, file = paste0(fasta_dir, "/", fasta_file_name), append=TRUE)
          setTxtProgressBar(pb1b, j)
        }
        
        if (nuc_search$count %% 250 != 0) {
          j <<- nuc_search$count - dif + 1
          tryCatch(write_me <- entrez_fetch(db = "nucleotide", web_history = cookie, rettype = "fasta", retstart = j, query_key = qk, retmax = dif), error = function(e) {pull2_again()})
          write(write_me, file = paste0(fasta_dir, "/", fasta_file_name), append=TRUE)
          setTxtProgressBar(pb1b, j)
        }
        close(pb1b)
      }
    }
    
  }
  rm(pull1, pull1_again, pull2, pull2_again, j)
}

library(rentrez)

EntrezDownload(txid = 9721, fasta_dir = fd, minimumSeq = 350, maximumSeq = NULL)

#### Run BUSCO ####

RunBusco <- function(fasta_dir, toast_dir, path_to_run_busco.py, path_to_orthoDB, threads = 1){
  reset_wd <- getwd()
  if (length(fasta_dir) == 0) {
    stop("We could not find any fasta files, exiting now\n")
  }
  
  species_fasta <- dir(fasta_dir)[grep(".fasta", dir(fasta_dir))]
  species_list <- gsub(".fasta", "", species_fasta)
  
  #check if any species have already been analysed with busco
  busco_dir <- paste0(toast_dir, "/busco_results/")
  busco_folders <- dir(busco_dir)[grep("run_", dir(busco_dir))]
  busco_species_prior <- gsub("run_", "", busco_folders)
  
  if (length(setdiff(species_list, busco_species_prior)) > 0) {
    differences <- setdiff(species_list, busco_species_prior)
    
    cat("\nPerforming BUSCO on", length(differences), "fasta file(s)\n\nSTEP 2 - RUNNING BUSCO\n")
    setwd(paste0(toast_dir, "/busco_results/"))
    pb2 <- txtProgressBar(min = 0, max = length(differences), style = 3)
    for (i in 1:length(differences)) {
      busco_fasta <- paste0(differences[i], ".fasta")
      system(paste0("python3 ", path_to_run_busco.py, " -c ", threads, " -i ", toast_dir,
                    "/fasta/", busco_fasta, " -o ", differences[i], " -l ", path_to_orthoDB,
                    " -m tran", " >> ", toast_dir, "busco_run_notes.txt"))
      setTxtProgressBar(pb2, i)
    }
    close(pb2)
  } else { cat("\nIt appears BUSCO has already been run on all given fasta files\n") }
  setwd(reset_wd)
}

RunBusco(fasta_dir = fd, toast_dir = td, path_to_run_busco.py = bs, path_to_orthoDB = od, threads = cpu)

#### Parse full busco tables ####

ParseBuscoResults <- function(busco_dir, ortho_dir, length_threshold = 0){
  
  if(length_threshold > 1) {
    stop("length_threshold should be a value between 0 and 1. A value of 0 will not set a required length for each busco id. A value of 1 will not include sequences that are less than the length of the busco id.")
  }
  
  busco_folders <- dir(busco_dir)[grep("run_", dir(busco_dir))] #refresh this here as new runs may have been created at STEP 2
  final_columns <- 'busco_id' #initialize the final column names
  
  busco_lengths <- read.table(paste0(ortho_dir, "/lengths_cutoff")) #lengths of each busco_id
  
  pb3 <- txtProgressBar(min = 0, max = length(busco_folders), style = 3)
  for (n in 1:length(busco_folders)){
    
    sub_busco <- gsub("run_", "", busco_folders[n])
    table_to_read <- paste(busco_dir, "/", busco_folders[n], "/full_table_", sub_busco, ".tsv", sep = "") #need to get species list in this function
    MyData <- read.delim(table_to_read, skip = 4, header = TRUE, sep = "\t", fill = TRUE, na.strings = "NA")
    
    if (exists("row_names") == FALSE){ #if it doesn't exist, create it
      row_names <- unique(MyData[,1]) #get the row names only the first pass through
    }
    
    if (exists("parsed_busco_df") == FALSE){ #if it doesn't exist, create it
      MyInitial <- data.frame(row_names, NA) #add the row_names to the final data
      matrix_holder <- MyInitial[,1]
      parsed_busco_df <- data.frame(matrix_holder)
    }
    
    final_columns <- append(final_columns, sub_busco)
    MyCleanData <- data.frame(NULL)
    targets <- unique(MyData[,1])
    
    for (i in 1:length(targets)){
      specific_species <- MyData[which(MyData[,1]==targets[i]),] #Shows sequences for each busco hit
      target_value <- max(specific_species$Score) #Select sequence with the better score
      row_needed <- specific_species[which(specific_species$Score==target_value),]
      
      if(!is.na(row_needed[1,1])) {
        current_length_threshold <- busco_lengths[which(busco_lengths[,1] == targets[i]), 4] * length_threshold #Create minimum length sequence should be to include
        
        if(row_needed$Length[1] < current_length_threshold) {
          row_needed[1,] <- NA
        }
      }
      
      MyCleanData <- rbind(MyCleanData,row_needed[1,]) #if scores above were identical, this grabs just the first one
      
    }
    setTxtProgressBar(pb3, n)
    parsed_busco_df <- cbind(parsed_busco_df, MyCleanData[,3]) #tosses out warnings when all are missing but this is fine to ignore
  }
  
  close(pb3)
  names(parsed_busco_df) <- final_columns #rename the columns based on the order they were parsed
  return(parsed_busco_df)
}

parsed_busco_results <- ParseBuscoResults(busco_dir = bd)

write.table(parsed_busco_results, file = paste0(td, "/parsed_busco_results.tsv"), sep = "\t", row.names = FALSE) #write it to a tsv inside the working directory

#### Extract BUSCO sequences ####

library(seqinr)

ExtractBuscoSeqs <- function(busco_table, fasta_dir, extract_dir){
  
  pb4 <- txtProgressBar(min = 0, max = length(names(busco_table)), style = 3)
  for (i in 2:length(busco_table[1,])) { #first column is busco.id so start the count at 2
    fasta_file <- read.fasta(file = paste0(fasta_dir, "/", colnames(busco_table)[i], ".fasta"), as.string = TRUE)
    setTxtProgressBar(pb4, i)
    for (j in 1:length(busco_table[,1])) {
      cell_value <- toString(busco_table[j,i])
      
      if (cell_value != "NA") {
        seq_to_write <- gsub("(.{60})", "\\1\n", fasta_file[[cell_value]][[1]])
        cat(">", colnames(busco_table)[i], "\n", seq_to_write, "\n", file = paste0(extract_dir, "/", busco_table[j,1], ".fasta"), append = TRUE, sep = "")
      }
    }
  }
  close(pb4)
}

ExtractBuscoSeqs(busco_table = parsed_busco_results, fasta_dir = fd, extract_dir = ed) #parsed_busco_results from previous step

#### Align using mafft ####

MafftOrientAlign <- function(extract_dir, mafft_dir, threads = 1, omit_invariants = T) {
  busco_id_matched_list <- list.files(path = extract_dir) #what happens if mafft tries to align a single sequence?
  pb5 <- txtProgressBar(min = 0, max = length(busco_id_matched_list), style = 3)
  for (i in 1:length(busco_id_matched_list)) {
    setTxtProgressBar(pb5, i)
    file_to_use <- busco_id_matched_list[i]
    system(paste0("mafft --adjustdirection --thread ", threads, " --quiet ", extract_dir, "/", file_to_use,
                  " > ", mafft_dir, "/", file_to_use))
    x <- readLines(paste0(mafft_dir, "/", file_to_use))
    y <- gsub("_R_", "", x)
    cat(y, file = paste0(mafft_dir, "/", file_to_use), sep = "\n")
    
    
    #if an aligned file that was just created is invariant, then remove that file
    if(omit_invariants == T) {
      num_fasta <- length(names(read.fasta(paste0(mafft_dir, "/", busco_id_matched_list[i]), as.string = T)))
      all_seqs <- c()
      
      for(j in 1:num_fasta) {
        all_seqs <- c(all_seqs, read.fasta(paste0(mafft_dir, "/", busco_id_matched_list[i]), as.string = T)[[j]][[1]])
      }
      
      if(length(unique(all_seqs)) == 1) {
        file.remove(paste0(mafft_dir, "/", busco_id_matched_list[i]))
      }
    }
  }
  close(pb5)
}

MafftOrientAlign(extract_dir = ed, mafft_dir = md, threads = cpu) #important as some sequences may be 3'->5' direction

#### Missing data check ####

mdf <- MissingDataTable(aligned_dir = ad)

#### Look at how much data is missing ####

colSums(is.na(mdf))

#### Remove species/samples that are missing too much data ####

ThresholdDataTable <- function(missing_df, threshold) {
  
  appended <- NULL
  max <- nrow(missing_df)
  summary <- colSums(is.na(missing_df)) #make a list of how many null variables are found
  for (i in 1:length(summary)){
    if (summary[i] < (max - threshold)){
      appended <- append(appended, summary[i])
    }
  }
  
  wanted <- names(appended)
  
  for (j in 1:length(wanted)){
    if (exists("threshold_df") == FALSE){ #if it doesn't exist, create it
      threshold_df <- missing_df[wanted[j]]
    } else {
      threshold_df <- cbind(threshold_df, missing_df[wanted[j]])
    }
  }
  return(threshold_df)
}

threshold_df <- ThresholdDataTable(missing_df = mdf, threshold = 1000)

# threshold_df can now serve as parsed_busco_table from previous steps to generate alignment #

ThresholdExtract <- function(aligned_dir, missing_df, threshold_fasta_folder){
  pb4 <- txtProgressBar(min = 0, max = nrow(missing_df), style = 3)
  for(j in 1:nrow(missing_df)){ #parse through each row, by column
    setTxtProgressBar(pb4, j) #change progress bar after each row (species)
    for(i in 1:ncol(missing_df)){
      cell_value <- toString(missing_df[j,i])
      if(cell_value != "NA" ){
        partial_fasta_name <- rownames(missing_df)[j]
        full_fasta_name <- paste0(partial_fasta_name, ".fasta")
        fasta_file <- read.fasta(file = paste0(aligned_dir, "/", full_fasta_name), as.string = TRUE)
        species_to_get <- colnames(missing_df)[i]
        seq_to_write <-  gsub("(.{60})", "\\1\n", fasta_file[[species_to_get]][[1]])
        cat(">", colnames(missing_df)[i], "\n", seq_to_write, "\n", file = paste0(threshold_fasta_folder, "/", full_fasta_name), append = TRUE, sep = "")
      }
    }
  }
  close(pb4)
}

ThresholdExtract(aligned_dir = ad, missing_df = threshold_df, threshold_fasta_folder = "/home/dustin/projects/holo_toast/threshold100")

#### Realign ####

MafftOrientAlign(extract_dir = "/home/dustin/projects/holo_toast/threshold1000", mafft_dir = "/home/dustin/projects/holo_toast/threshold1000/mafft_aligned", threads = cpu)

#### Table Partition ####

PartitionTable <- function(aligned_dir, missing_df){
  old_max <- 1
  cat("#nexus\nbegin sets;\n", file = "table.partition", append = TRUE)
  for (i in 1:length(row.names(missing_df))){
    if (max(missing_df[i,], na.rm = TRUE) > 0) {
      new_max <- old_max + max(missing_df[i,], na.rm = TRUE) #ignores any lines were all are NA
      new_start <- new_max - 1
      cat("charset ", row.names(missing_df)[i], " = ", old_max, "-", new_start, ";\n", sep = "", file = "table.partition", append = TRUE)
      old_max <- new_max
    }
  }
  cat("end;\n", file = "table.partition", append = TRUE)
}

PartitionTable(aligned_dir = "/home/dustin/projects/holo_toast/threshold1000/mafft_aligned", missing_df = threshold_df) #writes table.partition to working directory

#### Multi-gene alignment ####

SuperAlign <- function(aligned_dir, missing_df){
  num_species <- length(colnames(missing_df))
  old_max <- 0
  for (i in 1:length(row.names(missing_df))){
    if (max(missing_df[i,], na.rm = TRUE) > 0) {
      new_max <- old_max + max(missing_df[i,], na.rm = TRUE) #ignores any lines were all are NA
      old_max <- new_max
    }
  }
  
  max_align_length <- old_max
  header <- paste0(num_species, "\t", max_align_length)
  cat(file = "superalign.txt", header, "\n", append = TRUE, sep = '')
  
  for (i in 1:length(colnames(missing_df))){
    species_to_find <- colnames(missing_df)[i] #change to i later
    thing_to_append <- paste0(species_to_find, "\t")
    for (j in 1:nrow(missing_df)){
      if (is.na(missing_df[j,i]) == FALSE){ #grab the sequence
        fasta_file <- read.fasta(paste0(aligned_dir, "/", row.names(missing_df[j,]), ".fasta"), as.string = TRUE)
        target_species <- match(species_to_find, names(fasta_file))
        seq_to_append <- fasta_file[[target_species]][[1]]
        thing_to_append <- paste0(thing_to_append, seq_to_append)
      }
      if (is.na(missing_df[j,i]) == TRUE){
        dashes <- max(missing_df[j,], na.rm = TRUE) #add in this number of dashes
        dashes_to_append <- strrep("-", dashes)
        thing_to_append <- paste0(thing_to_append, dashes_to_append)
      }
    }
    cat(file = "superalign.txt", thing_to_append, "\n", append = TRUE, sep = '')
  }
}

SuperAlign(aligned_dir, missing_df = mdf)

#############################
#############################
#### VISUALIZE YOUR DATA ####
#############################
#############################

### have a look at the biological processes assigned to the orthoDB BUSCO Ids via Gene Ontology

library(GSEABase)

ToastGoSlim <- function(orthogroup_info, obo, perspective = "BP"){
  orthogroup_info <- read.delim(orthogroup_info, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
  pbcounter <- 0 #will help keep track of progression in progress bar through loops
  pb7 <- txtProgressBar(min = 0, max = nrow(orthogroup_info), style = 3)
  for (i in 1:nrow(orthogroup_info)){
    go_ids <- unlist(strsplit(as.character(orthogroup_info[i,"BiologicalProcesses"]), ";"))#gene ontology IDs imported from orthoDB info file
    if (length(go_ids) == 1 && is.na(go_ids) == TRUE){ #NA is still length 1; throws a bunch of warnings if it has to evaluate chr of length > 1
      if (exists("empty") == FALSE){ #start the variable if it doesn't exist
        empty <- unlist(as.character(orthogroup_info[i, "OrthoGroupID"]))
      } else { #make a list of empty orthogroups to append to the end
        empty <- append(empty, unlist(as.character(orthogroup_info[i, "OrthoGroupID"])))
      }
      pbcounter <- pbcounter + 1 #keep track of progress through the three perspectives loops
      next()
    } else {
      myCollection <- GOCollection(go_ids)
      slim <- getOBOCollection(obo)
      slimmed <- goSlim(myCollection, slim, perspective)
      row.names(slimmed) <- slimmed$Term
      if (exists("headers") == FALSE) { #start the variable if it doesn't exist
        headers <- unlist(as.character(orthogroup_info[i, "OrthoGroupID"]))
      } else {
        headers <- append(headers, unlist(as.character(orthogroup_info[i, "OrthoGroupID"])))
      }
      if (exists("appended") ==  FALSE) { #start the variable if it doesn't exist
        appended <- slimmed[,1]
      } else {
        appended <- cbind(appended, slimmed[,1])
      }
      pbcounter <- pbcounter + 1 #keep track of progress through the three perspectives loops
    }
    setTxtProgressBar(pb7, pbcounter)
  }
  colnames(appended) <- headers
  #need to merge back in the orthogroup_info[i,"BiologicalProcesses"] that were empty
  zeros <- nrow(appended) * length(empty) #will help populate the emptry matrix
  empty_matrix <- matrix(rep(0, zeros), nrow = nrow(appended)) #create an matrix full of zeros with
  colnames(empty_matrix) <- empty                              #colnames to merge with appended
  merger <- cbind(appended, empty_matrix) #this is the desired matrix
  row.names(merger) <- slimmed$Term
  return(merger)
  close(pb7)
}

BP <- ToastGoSlim(orthogroup_info = ortho, obo = go_oro, perspective = "BP") #converts BUSCO_IDs into a dataframe of Biological Process GoSlim Terms
#MF and CC are currently not assigned in orthoDB and therefore not supported in TOAST
#Visualize the overlapping GoSlim terms using the nifty package UpSetR
# check out their project here: https://github.com/hms-dbmi/UpSetR
if ("UpSetR" %in% rownames(installed.packages()) == FALSE) {
  install.packages("UpSetR")
}

library(UpSetR)
#Create intersection plot
#format BP data set to plot using the upset() function
rownames(BP)[21] <- "carb. derivative metabolic" #original name too long
to_plot <- t(BP) #transpose the data
to_plot[to_plot > 0] <- 1 #remove comparisons with no overlap
to_plot <- data.frame(to_plot) # convert to dataframe
to_plot <- cbind(rownames(to_plot), to_plot) #set rownames to be the first row, remove rownames later
colnames(to_plot)[1] <- "Identifier"
rownames(to_plot) <- NULL #remove rownames

#format BP data set to plot
to_plot[,2:ncol(to_plot)] <- sapply(to_plot[2:ncol(to_plot)], as.integer)

#create intersection plot
upset(to_plot, order.by = "freq", sets = colnames(to_plot)[2:ncol(to_plot)])

#### Create an occupancy matrix plot to display which Busco Ids were found in each sample

VisualizePBR <- function(busco_dir){
  
  busco_folders <- dir(busco_dir)[grep("run_", dir(busco_dir))] #refresh this here as new runs may have been created at STEP 2
  final_columns <- 'busco_id' #initialize the final column names
  pb3 <- txtProgressBar(min = 0, max = length(busco_folders), style = 3)
  for (n in 1:length(busco_folders)){
    sub_busco <- gsub("run_", "", busco_folders[n])
    table_to_read <- paste(busco_dir, "/", busco_folders[n], "/full_table_", sub_busco, ".tsv", sep = "") #need to get species list in this function
    MyData <- read.delim(table_to_read, skip = 4, header = TRUE, sep = "\t", fill = TRUE, na.strings = "NA")
    
    if (exists("row_names") == FALSE){ #if it doesn't exist, create it
      row_names <- unique(MyData[,1]) #get the row names only the first pass through
    }
    
    if (exists("parsed_busco_df") == FALSE){ #if it doesn't exist, create it
      MyInitial <- data.frame(row_names, NA) #add the row_names to the final data
      matrix_holder <- MyInitial[,1]
      parsed_busco_df <- data.frame(matrix_holder)
    }
    
    final_columns <- append(final_columns, sub_busco)
    MyCleanData <- data.frame(NULL)
    targets<-unique(MyData[,1])
    
    for (i in 1:length(targets)){
      specific_species <- MyData[which(MyData[,1]==targets[i]),]
      target_value <- max(specific_species$Score)
      row_needed <- specific_species[which(specific_species$Score==target_value),]
      MyCleanData <- rbind(MyCleanData,row_needed[1,]) #if scores above were identical, this grabs just the first one
    }
    setTxtProgressBar(pb3, n)
    parsed_busco_df <- cbind(parsed_busco_df,MyCleanData[,2]) #[,2] should grab the status (missing, fragmented, complete, duplicated)
  }
  
  close(pb3)
  names(parsed_busco_df) <- final_columns #rename the columns based on the order they were parsed
  
  parsed_busco_df[] <- lapply(parsed_busco_df, as.character)
  #change the values in the dataframe to numeric
  #        which(parsed_busco_df == "Fragmented") <- 0.5
  #        which(parsed_busco_df == "Complete") <- 1
  #        which(parsed_busco_df == 'Duplicated') <- 1
  #        which(parsed_busco_df == "NA") <- 0
  return(parsed_busco_df)
}

VPBR <- VisualizePBR(busco_dir = bd) #creates a data.frame, similar to the function ParseBuscoResults()

dim(VPBR) #check out the dimensions

vpbr <- as.matrix(VPBR[,2:length(colnames(VPBR))]) #convert VPBR data.frame to vpbr matrix
vpbr <- gsub("Fragmented", 0.5, vpbr)              # format the data frame into a matrix and replace Fragmented,
vpbr <- gsub("Complete", 1.0, vpbr)                # Complete, Duplicated and NA (ie missing) with numeric values
vpbr <- gsub("Duplicated", 1.0, vpbr)
vpbr[is.na(vpbr)] <- 0

class(vpbr) <- "numeric" #image function below needs numeric or logical values to plot

vpbr <- vpbr[,c("busco_Cambac_brain.fasta", "busco_Camdro_brain.fasta",
                "busco_Cambac_hypothalamus.fasta", "busco_Camdro_hypothalamus.fasta",
                "busco_Cambac_kidney.fasta", "busco_Camdro_kidney.fasta",
                "busco_Cambac_liver.fasta", "busco_Camdro_liver.fasta",
                "busco_Cambac_lung.fasta", "busco_Camdro_lung.fasta",
                "busco_Cambac_muscle.fasta", "busco_Camdro_muscle.fasta",
                "busco_Cambac_skin.fasta", "busco_Camdro_skin.fasta",
                "busco_Cambac_testis.fasta", "busco_Camdro_testis.fasta")]

# Create y axis labels
yLabels <- gsub("busco_", "", colnames(vpbr))#remove the leading busco_ and trailing .fasta
yLabels <- gsub(".fasta", "", yLabels)

ColorRamp <- gray.colors(3, start = 1, end = 0, gamma = 2.2, rev = FALSE) #white = absent, gray = fragmented, black = duplicated/complete
Colors <- c("#FFFFFF", "#BABABA", "#57A0D3")
# Set layout
par(mar = c(5,15,2.5,1), font = 2) #mar=C(bot, left, top, righ 45-t)

# Make the plot
image(1:nrow(vpbr), 1:ncol(vpbr), vpbr,
      col=Colors, xlab="", ylab="",
      axes=FALSE, main= NA)

# Annotate the y axis
box()
axis(side = 2, at=seq(1,length(yLabels),1), labels=yLabels, las= 1,
     cex.axis=1)

### Create Stacked Bar Graph showing how many orthologs were found belonging to the 21 GoTerms

# create color palette:
library(RColorBrewer)
coul <- brewer.pal(8, "Pastel1")
coul <- c(coul, brewer.pal(8, "Pastel2"))

data_percentage <- as.matrix(sorted_pie_df/rowSums(sorted_pie_df)*100)

class(data_percentage)

# Make a stacked barplot--> it will be in %!
end_point = 0.5 + nrow(sorted_pie_df) + nrow(sorted_pie_df)-1
par(mar = c(15,5,1,15))
barplot(t(data_percentage), col=coul , border="white", space = 1.0, xlab = "", ylab = "Percent",
        main = "Percent Sample within GoTerm Category",
        axisnames = F)
text(seq(1.5,end_point,by=2), par("usr")[3]-0.25,
     srt = 60, adj= 1, xpd = TRUE,
     labels = paste(rownames(sorted_pie_df)), cex=0.9)

legend(x = 45, y = 100, legend = rev(colnames(sorted_pie_df)), fill = rev(coul), xpd = 1)

