#!/bin/env Rscript

sampleSize <- 10000
count <- 5
p <- 0.01

for(i in 1:sampleSize) {
    v <- rgeom(count,p)
    s <- sum(v)
    cat(s,"\n")
    #sample[i] <- s
}


