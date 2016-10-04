#!/bin/env Rscript

args <- commandArgs(TRUE)
if(length(args)!=2) {
  cat("usage: P N\n");
  q(status=1)
}
p <- as.numeric(args[1])
count <- as.numeric(args[2])
sampleSize <- 10000

for(i in 1:sampleSize) {
    v <- rgeom(count,p)
    s <- sum(v)
    cat(s,"\n")
    #sample[i] <- s
}


