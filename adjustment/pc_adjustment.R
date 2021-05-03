
library(sva)
library(WGCNA)
library(data.table)
library(tidyverse)
library(variancePartition)
library('doParallel')
library(readxl)
library(plyr)
library(Hmisc)
library(confeti)
library(corrplot)
library(RColorBrewer)
library(gplots)
library(biomaRt)
library(factoextra)
library(gsubfn)

expr <- args[1]

mod=matrix(1,nrow=dim(expr.t)[1],ncol=1)
colnames(mod)="Intercept"
nsv=num.sv(t(expr.t),mod, method = "be")

print(paste("Number of PCs estimated to be removed:", nsv))

# PC residualization of GE matrix
exprs_corrected = sva_network(as.matrix(expr.t), nsv)
dat <- as.data.frame(exprs_corrected) %>% rownames_to_column("SampleID")
