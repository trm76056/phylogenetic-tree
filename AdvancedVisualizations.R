#this script will run ToastGoSlim to identify GO terms for the data, generate the UpSetR plot for visualizing that,

#' Assign Go_Slim IDs to BUSCO IDs
#'
#' this function takes the Gene Ontology information contained in the orthoDB and generalizes it
#'         into as as few GO terms as possible given three perspectives: Molecular function,
#'         Biological Process and Cellular Component, each returned as a dataframe
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @keywords toast missing transcript sequence DNA phylogeny fasta busco ortholog gene ontology go slim
#' @param orthogroup_info assing the uncompressed ontology file found in orthoDB/info/*orthogroup_info.txt.gz
#' @param obo download desired ontology information (.obo file) from http://geneontology.org/docs/download-ontology/#go_obo_and_owl
#'            these are consistently updated so make sure to grab the newest
#'            we recommend the GO slim AGR subset (goslim_agr.obo) which can be obtained using the command
#'            wget http://current.geneontology.org/ontology/subsets/goslim_agr.obo
#' @param perspective options are BP(biological process), MF(Molecular function), or CC(cellular component)
#' @import GOstats GSEABase BiocManager
#' @export
#' @examples
#' BP <- ToastGoSlim(orthogroup_info = "path/to/orthoDB/info/*orthogroup_info.txt.gz", obo = "path/to/goslim_agr.obo", perspective = "BP")

library("GOstats")

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

#go_oro = "/scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/goslim_agr.obo"
#ortho = "/scratch/trm76056/PhyloTreeData/My_TOAST_Project_odb9/metazoa_odb9/info/metazoa_33208_OrthoDB9_orthogroup_info.txt.gz"
go_oro = "~/Documents/Research/My_TOAST_Project/goslim_agr.obo"
ortho = "~/Downloads/metazoa_odb9/info/metazoa_33208_OrthoDB9_orthogroup_info.txt"

BP <- ToastGoSlim(orthogroup_info = ortho, obo = go_oro, perspective = "BP")

if ("UpSetR" %in% rownames(installed.packages()) == FALSE) {
 install.packages("UpSetR")
}

library(UpSetR)
#Create intersection plot
#format BP data set to plot using the upset() function
#rownames(BP)[21] <- "carb. derivative metabolic" #original name too long
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

#create occupancy matrix plot to disply which busco ids were found in each sample

#' Explore Global Patterns of Missing Data
#'
#' Returns a dataframe that can be used to generated graphics of missing data patterns
#' @author Dustin J Wcisel, \email{djwcisel@@ncsu.edu}
#' @author James Thomas Howard, \email{jthowar3@@ncsu.edu}
#' @author Jeffrey A Yoder, \email{jayoder@@ncsu.edu}
#' @author Alex Dornburg, \email{dornburgalex@@gmail.com}
#' @param tsv A data frame of missing data coverage such as the generated “parsed_busco_results.tsv”
#' @param threshold Cutoff value for minimum number of sequences allowed per taxon
#' @keywords toast missing decisiveness sequence DNA phylogenomics
#' @export
#' @examples
#' VPBR <- isualizePBR(busco_dir = bd)

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

bd <- "~/Documents/Research/My_TOAST_Project/busco_results/"

VPBR <- VisualizePBR(busco_dir = bd) #creates a data.frame, similar to the function ParseBuscoResults()
dim(VPBR) #check out the dimensions
vpbr <- as.matrix(VPBR[,2:length(colnames(VPBR))]) #convert VPBR data.frame to vpbr matrix
vpbr <- gsub("Fragmented", 0.5, vpbr) # format the data frame into a matrix and replace Fragmented,
vpbr <- gsub("Complete", 1.0, vpbr) # Complete, Duplicated and NA (ie missing) with numeric values
vpbr <- gsub("Duplicated", 1.0, vpbr)
vpbr[is.na(vpbr)] <- 0
class(vpbr) <- "numeric" #image function below needs numeric or logical values to plot

vpbr<- vpbr[,c("busco_acutilabella.ORP","busco_albomicans.ORP","busco_angularis.ORP","busco_arawakana.ORP","busco_bizonata.ORP","busco_brachynephros.ORP","busco_cardini.ORP","busco_cardinoides.ORP","busco_deflecta.ORP","busco_dunni.ORP","busco_falleni.ORP","busco_funebris.ORP","busco_guarani.ORP","busco_guttifera.ORP","busco_immigrans.ORP","busco_macroptera.ORP","busco_macrospina.ORP","busco_multispina.ORP","busco_munda.ORP","busco_neotestacea.ORP","busco_nigrodunni.ORP","busco_nigromaculata.ORP","busco_occidentalis.ORP","busco_orientacea.ORP","busco_pallidipennis.ORP","busco_palustris.ORP","busco_parthenogenetica.ORP","busco_phalerata.ORP","busco_polymorpha.ORP","busco_putrida.ORP","busco_quinaria.ORP","busco_recens.ORP","busco_similis.ORP","busco_subbadia.ORP","busco_suboccidentalis.ORP","busco_subpalustris.ORP","busco_subquinaria_Co.ORP","busco_subquinaria_In.ORP","busco_sulfurigaster.ORP","busco_tenebrosa.ORP","busco_testacea.ORP","busco_transversa.ORP","busco_tripunctata.ORP")]

# Create y axis labels
yLabels <- gsub("busco_", "", colnames(vpbr))#remove the leading busco_ and trailing .fasta
yLabels <- gsub(".fasta", "", yLabels)
ColorRamp <- gray.colors(3, start = 1, end = 0, gamma = 2.2, rev = FALSE) #white =absent, gray = fragmented, black = duplicated/complete
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

#visualize the percent of each sample found within a GOterm category

# create color palette:
library(RColorBrewer)
coul <- brewer.pal(8, "Pastel1")
coul <- c(coul, brewer.pal(8, "Pastel2"))
data_percentage <- as.matrix(sorted_pie_df/rowSums(sorted_pie_df)*100)
class(data_percentage)
# Make a stacked barplot--> it will be in %!
end_point = 0.5 + nrow(sorted_pie_df) + nrow(sorted_pie_df)-1
par(mar = c(15,5,1,15))
barplot(t(data_percentage), col=coul , border="white", space = 1.0, xlab = "", ylab =
"Percent",
 main = "Percent Sample within GoTerm Category",
 axisnames = F)
text(seq(1.5,end_point,by=2), par("usr")[3]-0.25,
 srt = 60, adj= 1, xpd = TRUE,cex=0.9,
 labels = paste(rownames(sorted_pie_df)))
legend(x = 45, y = 100, legend = rev(colnames(sorted_pie_df)), fill = rev(coul), xpd =
1)
