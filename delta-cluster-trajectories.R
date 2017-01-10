#!/data/reddylab/software/miniconda2/bin/Rscript
options(width=500)
data <- read.table("trajectories.txt")
elements <- data[1]
data <- data[2:12]
trinary <- ifelse(data>1000,1,ifelse(data< -1000,-1,0))
clusters <- kmeans(trinary,10,iter.max=100)
elements <- cbind(elements,clusters$cluster)
write.table(elements,file="cluster-membership.txt",sep="\t", col.names = F, row.names = F)
write.table(clusters$centers,file="clusters.txt",sep="\t", col.names = F, row.names = F)
write.table(clusters$size,file="cluster-sizes.txt",sep="\t", col.names = F, row.names = F)


