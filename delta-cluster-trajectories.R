#!/data/reddylab/software/miniconda2/bin/Rscript
options(width=500)

args <-commandArgs(TRUE)
if(length(args)!=4) {
    cat("usage: trajectories.txt membership.txt clusters.txt sizes.txt\n")
    q(status=1)
}
infile <- args[1]
membershipFile <- args[2]
clustersFile <- args[3]
sizesFile <- args[4]

#data <- read.table("trajectories-onepath.txt")
data <- read.table(infile)
elements <- data[1]
data <- data[2:12]
trinary <- ifelse(data>1000,1,ifelse(data< -1000,-1,0))
clusters <- kmeans(trinary,10,iter.max=100)
elements <- cbind(elements,clusters$cluster)
write.table(elements,file=membershipFile,sep="\t", col.names = F, row.names = F)
write.table(clusters$centers,file=clustersFile,sep="\t", col.names = F, row.names = F)
write.table(clusters$size,file=sizesFile,sep="\t", col.names = F, row.names = F)

#write.table(elements,file="cluster-membership.txt",sep="\t", col.names = F, row.names = F)
#write.table(clusters$centers,file="clusters.txt",sep="\t", col.names = F, row.names = F)
#write.table(clusters$size,file="cluster-sizes.txt",sep="\t", col.names = F, row.names = F)


