#!/data/reddylab/software/miniconda2/bin/Rscript
options(width=500)

data <- read.table("feature-subsets.txt",header=F,sep=" ",dec=".",strip.white=TRUE);
plotTop <- 0.9
pdf("barplot.pdf")
par(mfrow=c(2,1))
barCenters <- barplot(data$V2,width=1,axes = TRUE,names.arg = data$V1,ylim = c(0, plotTop),xlim=c(0,22),col="red")
d1 <- data$V2 - data$V3/5 # * 2
d2 <- data$V2 + data$V3/5 # * 2
segments(barCenters, d1, barCenters, d2, lwd = 1.5)
arrows(barCenters, d1, barCenters, d2, lwd = 1.5, angle = 90, code = 3, length = 0.05)
dev.off()
