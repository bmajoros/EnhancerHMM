#!/data/reddylab/software/miniconda2/bin/Rscript
# #!/bin/env Rscript
args <-commandArgs(TRUE)
if(length(args)!=2) {
    cat("usage: edgeR-intput.txt out.parms\n")
    q(status=1)
}

library(edgeR)
mat <- args[1]
y_mat_and_dispersion_out <- args[2]

mat <- read.table(mat, header = T, row.names=1, sep='\t')
factors <- factor(c("DNA","ETH","ETH","ETH","DEX","DEX","DEX"))

model_design <- model.matrix(~factors)

print(paste('Observations, =', min(dim(model_design))))
print(paste('Parameters to estimate, =',qr(model_design)$rank))

####
y <- DGEList(counts = mat, group = factors)

y <- calcNormFactors(y, design=model_design)
y <- estimateGLMCommonDisp(y, design=model_design)
y <- estimateGLMTrendedDisp(y, design=model_design) 
y <- estimateGLMTagwiseDisp(y, design=model_design)

saveRDS(y, y_mat_and_dispersion_out)


