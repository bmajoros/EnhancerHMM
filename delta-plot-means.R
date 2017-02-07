#!/data/reddylab/software/miniconda2/bin/Rscript
options(width=500)

args <-commandArgs(TRUE)
if(length(args)!=2) {
    cat("usage: in.means out.pdf\n")
    q(status=1)
}
infile <- args[1]
outfile <- args[2]

maxY <- 2.5
minY <- -1.4
#maxY <- 2.0
#minY <- -2.5
#maxY <- 1.5
#minY <- -1.5
pointtype <- 16
linewidth <- 2

#df0 <- read.table("means-t0.txt")
#df0 <- as.data.frame(t(df0))
#df3 <- read.table("means-t3.txt")
#df3 <- as.data.frame(t(df3))
#df <- as.data.frame(as.matrix(df3)-as.matrix(df0))
#dfmax <- max(df)
#dfmin <- min(df)

df <- read.table(infile)
df <- as.data.frame(t(df))
dfmax <- maxY
dfmin <- minY
#cat(paste(dfmin,"\t",dfmax))
df <- cbind(df,c(1,2,3,4,5))

pdf(outfile)
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


