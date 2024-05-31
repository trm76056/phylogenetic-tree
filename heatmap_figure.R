#This script will be used to help create the heatmap figure for the sensitivity analysis

setwd("~/Documents/Research/My_TOAST_Project/")

#load required packages

library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
#library(ggtree)

#set key parameters
#loci_cutoffs <- c("25","50","75","100","125","150","175","200","225","250","275",
                  #"300","325","350","375","400","425","450","475","500","525","550","575",
                  #"600","625","650","675","700","725","750","775","800","825","850","875",
                  #"900","925","950")
loci_cutoffs <- c("50","100","150","200","250","300",
                  "350","400","450","500","550","600",
                  "650","700","750","800","850","900","950")
bins_levels <- c("2","3","4","5","6","7","8")
heatmap_colors <- c("white","khaki1", "cyan","mediumblue","black")
heatmap_bins <- c(0, 25, 50, 75, 100)
#nodes <- c(1:47) #this is how many nodes you have (48 species = 47 nodes) for astral tree
nodes <- c(1:46) #46 for raxml tree idk why

#tree <- read.nexus("Astral_bootstrapoutput_all.nexus") #all trees must have same node labels, might need to reroot
tree <- read.nexus("RAxML_bootstrapoutput_all.nexus")

#get BS node values from each tree
dat <- lapply(tree, '[[', "node.label")
dat <- data.frame(dat)
dat <- tibble::as_tibble(dat)

#create new column names
#incl_3 <- c('3_61','3_122','3_183','3_244','3_305','3_366','3_427','3_489')
#incl_3 <- c('3_75','3_150','3_225','3_300','3_375','3_375','3_415', '3_490)
incl_3 <- c('2_50','2_100','2_200','2_250','2_300','2_350','2_450','2_500')
#incl_4 <- c('4_65','4_130','4_195','4_260','4_325','4_390','4_455','4_517')
#incl_4 <- c('4_75','4_125','4_200','4_250','4_325','4_400','4_450','4_525')
incl_4 <- c('3_50','3_150','3_200','3_250','3_300','3_400','3_450','3_500')
#incl_5 <- c('5_87','5_174','5_261','5_348','5_435','5_522','5_609','5_696')
#incl_5 <- c('5_100','5_175','5_250','5_350','5_425','5_525','5_600','5_700')
incl_5 <- c('4_100','4_150','4_250','4_350','4_450','4_500','4_600','4_700')
#incl_6 <- c('6_108','6_216','6_324','6_432','6_540','6_648','6_756','6_866')
#incl_6 <- c('6_100','6_225','6_325','6_425','6_550','6_650','6_750','6_875')
incl_6 <- c('5_100','5_200','5_300','5_400','5_550','5_650','5_750','5_850')
#incl_7 <- c('7_116','7_232','7_348','7_464','7_580','7_696','7_812','7_930')
#incl_7 <- c('7_125','7_225','7_350','7_450','7_575','7_700','7_825','7_925')
incl_7 <- c('6_100','6_250','6_350','6_450','6_600','6_700','6_800','6_950')
#incl_8 <- c('8_118','8_236','8_354','8_472','8_590','8_708','8_826','8_945')
#incl_8 <- c('8_825','8_225','8_350','8_475','8_600','8_700','8_825','8_950')
incl_8 <- c('7_100','7_250','7_350','7_450','7_600','7_700','7_800','7_950')
#incl_9 <- c('9_118','9_236','9_354','9_472','9_590','9_708','9_826','9_947')
#incl_9 <- c('9_125','9_225','9_350','9_475','9_600','9_700','9_825','9_950')
incl_9 <- c('8_100','8_250','8_350','8_450','8_600','8_700','8_800','8_950')
#incl_10 <- c('10_118','10_236','10_354','10_472','10_590','10_708','10_826','10_947')
#incl_10 <- c('10_125','10_225','10_350','10_475','10_600','10_700','10_825','10_950')
incl_10 <- c('9_100','9_250','9_350','9_450','9_600','9_700','9_800','9_950')
incl <- c(incl_3,incl_4,incl_5,incl_6,incl_7,incl_8,incl_9,incl_10)
colnames(dat) <- incl
dat <- data.frame(dat)

#reorganize data into proper format
dat <- gather(dat, key = "set", value = "bs")
bs <- dat$bs
bs <- replace(bs, bs=="", "100.0")
dat <- subset(dat, select = -c(bs))
dat <- cbind(dat,bs)
dat <- separate(dat, set, sep = "_", into = c("bins", "loci")) #splits NUMbin_NUMloci column into two based on underscore
dat <- cbind(nodes, dat) #add columns
dat <- dat %>% drop_na()
dat$bins<-gsub("X","",as.character(dat$bins))
#dat$bins <- gsub(pattern = "\\D", replacement = "", dat$bins) #gets rid of text in column
#dat$loci <- gsub(pattern = "\\D", replacement = "", dat$loci) #gets rid of text in column

#make figure with node values and bs values to double check the correct bs value is placed with correct node
#ggtree(trees) + geom_text2(aes(subset=!isTip, label=label)) #get BS values
#ggtree(trees) + geom_text2(aes(sebset=!isTip, label=node)) #get node labels

#generate heatmap for node
for (i in 1:length(nodes)) {
  #dat <- within(dat, loci <- factor(loci, levels = loci_cutoffs))
  #dat <- within(dat, bins <- factor(bins, levels = bins_levels))
  node <- dat[which(dat$nodes==i),] #make a figure for each individual node
  node <- data.frame(apply(node,2,as.numeric))
  #name <- paste('node', i, '.png', sep='')
  name <- paste('ASTRAL_node',i,'.png',sep='')
  #p <- ggplot(node, aes(loci,bins))+geom_tile(aes(fill = bs),width=150) + scale_fill_gradientn(limits=c(0,100),colours = heatmap_colors, values=heatmap_bins) + scale_x_continuous(limits=c(0,950))
  p <- ggplot(node, aes(loci,bins))+geom_tile(aes(fill = bs),width=150) + scale_fill_gradient(limits=c(0,100), low="white",high="black") + scale_x_continuous(limits=c(0,950)) +
      theme_minimal() + theme(plot.background = element_blank(), legend.position = 'none',
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.border = element_rect(fill = NA, colour = "black", size = 3),
                              axis.title.x = element_blank(), axis.title.y = element_blank(),
                              axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                              axis.text.y = element_blank(), axis.ticks.y = element_blank())
  png(file=name)
  #pdf(file = name)
  print(p)
  #ggsave(name, plot = 'p', device='png')
  dev.off()
}

#legend for sensitivity figure
i <- 13
node <- dat[which(dat$nodes==i),] #make a figure for each individual node
node <- data.frame(apply(node,2,as.numeric))
name <- paste('Sensitivity_Legend13',i,'.png',sep='')
p <- ggplot(node, aes(loci,bins))+geom_tile(aes(fill = bs),width=150) + scale_fill_gradient(limits=c(0,100), low="white",high="black") + scale_x_continuous(limits=c(0,950)) + theme_minimal() + theme(text = element_text(size = 40)) + labs(title="Legend",x="Number of Loci",y="Number of Bins")
png(file=name)
print(p)
dev.off()



