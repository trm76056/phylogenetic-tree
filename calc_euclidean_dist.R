#pairwise tree distance in MDS using R
#this is using the output from treeCMP

bins <- read.table(file="bins2345678910_treecmp.out", header=T)
gene_names <- read.table("bins2345678910_genenames.txt")
n <- max(table(bins$Tree1))+1
res <- lapply(with(bins, split(Triples, bins$Tree1)), function(x) c(rep(NA, n-length(x)),x))
res <- do.call("rbind", res)
res <- rbind(res, rep(NA,n))
res <- as.dist(t(res))

#set up mds
mds <- cmdscale(res, k=2)
plot(mds[,1], mds[,2], xlab="", ylab="", axes=TRUE,
     main="cmdscale (stats)")
#text(mds[,1], mds[,2], labels(res), cex = 0.9, xpd = TRUE) #use if you want tree names

#calculate euclidean distance from origin
#add origin to mds data
mds_origin <- rbind(c(0,0), mds)
#calculate euclidean distance of each point from the origin
bins_eucl <- as.matrix(dist(mds_origin, method="euclidean"))[,1]
#write the add, add the gene names and use that for rankings
bins_eucl <- data.frame(bins_eucl)[-1,]
bins_eucl <- data.frame(bins_eucl)
bins_eucl <- cbind(gene_names, bins_eucl)
#library(tidyverse)
bins_eucl_sort <- bins_eucl %>% arrange(bins_eucl)
write.csv(bins_eucl_sort, file="bins2345678910_euclidean_dist.csv")
