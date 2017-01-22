#!/data/reddylab/software/miniconda2/bin/Rscript
options(width=500)
#data <- read.table("trajectories.txt")
data <- read.table("trajectories-num-peaks.txt")
elements <- data[1]
data <- data[2:12]
trinary <- data
clusters <- kmeans(trinary,10,iter.max=100)
elements <- cbind(elements,clusters$cluster)
write.table(elements,file="numpeaks-cluster-membership.txt",sep="\t", col.names = F, row.names = F)
write.table(clusters$centers,file="numpeaks-clusters.txt",sep="\t", col.names = F, row.names = F)
write.table(clusters$size,file="numpeaks-cluster-sizes.txt",sep="\t", col.names = F, row.names = F)


