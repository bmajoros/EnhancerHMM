#!/data/reddylab/software/miniconda2/bin/Rscript
options(width=500)

pointtype <- 16
linewidth <- 2

df <- read.table("means.txt")
df <- as.data.frame(t(df))
dfmax <- max(df)
dfmin <- min(df)
#cat(paste(dfmin,"\t",dfmax))
df <- cbind(df,c(1,2,3,4,5))

pdf("means.pdf")
par(mfrow=c(2,1))

# DNase
plot(df[c(6,1)],ylim=c(dfmin,dfmax),pch=pointtype,col="red")
lines(df[c(6,1)],lwd=linewidth,col="red")

# P300
points(df[c(6,2)],pch=pointtype,col="orange")
lines(df[c(6,2)],lwd=linewidth,col="orange")

# H3K27ac
points(df[c(6,3)],pch=pointtype,col="blue")
lines(df[c(6,3)],lwd=linewidth,col="blue")

# H3K4me1
points(df[c(6,4)],pch=pointtype,col="cyan")
lines(df[c(6,4)],lwd=linewidth,col="cyan")

# H3K4me2
points(df[c(6,5)],pch=pointtype,col="navy")
lines(df[c(6,5)],lwd=linewidth,col="navy")

dev.off()


