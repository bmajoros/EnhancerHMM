#!/data/reddylab/software/miniconda2/bin/Rscript
# #!/bin/env Rscript
options(max.print = 99999999)
options(width=500)

args <-commandArgs(TRUE)
if(length(args)!=3) {
    cat("usage: edgeR-full-inputs.txt egdeR-thinned.txt out.parms\n")
    q(status=1)
}

library(edgeR)
fullmat <- args[1]
mat <- args[2]
y_mat_and_dispersion_out <- args[3]

mat <- read.table(mat, header = T, row.names=1, sep='\t')
fullmat <- read.table(fullmat,header=T,row.names=1,sep='\t')

# These need to be alphabetic because model.matrix() sorts them...
factors <- factor(c("A_DNA","B_ETH","B_ETH","B_ETH","C_DEX","C_DEX","C_DEX"))
model_design <- model.matrix(~factors)

# Estimate common dispersion from the thinned data
y <- DGEList(counts = mat, group = factors)
y <- calcNormFactors(y, design=model_design)
y <- estimateGLMCommonDisp(y, design=model_design)
y <- estimateGLMTrendedDisp(y, design=model_design)
common <- y$common.dispersion
trended <- y$trended.dispersion

# Up-sample the trended values to the full matrix
n <- nrow(fullmat)
fulltrend <- 1:n
for(i in 1:n) {
    fulltrend[[i]] <- trended[[as.integer(i/350)+1]]
}

# Estimate tagwise dispersion from the full data
y <- DGEList(counts = fullmat, group = factors)
y <- calcNormFactors(y, design=model_design)
y$common.dispersion <- common
y$trended.dispersion <- fulltrend
y <- estimateGLMTagwiseDisp(y, design=model_design)
saveRDS(y, y_mat_and_dispersion_out)

# Fit GLM model and perform likelihood ratio tests
fit <- glmFit(y, model_design)
cat("writing eth")
sink("edgeR-eth.txt")
de <- as.data.frame(topTags(glmLRT(fit, coef=2), n=dim(y)[1]))
print(de)
sink()
cat("writing dex")
sink("edgeR-dex.txt")
de <- as.data.frame(topTags(glmLRT(fit, coef=3), n=dim(y)[1]))
print(de)
sink()

