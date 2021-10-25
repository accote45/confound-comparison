


# Method = ICA FDR (ICA followed by FDR estimation)

# NOTE: original publication used the python implementation of this program
# detailed reference for ICA : http://www.cs.jhu.edu/~ayuille/courses/Stat161-261-Spring14/HyvO00-icatut.pdf

# source matrix: contains for every module the evidence that a certain gene belongs to the module
# mixing matrix: contains the individual contributions of a sample to every module

# X = pre-processed matrix
# K = pre-whitening matrix that projects data onto the first n.comp PCs
# W = estimated un-mixing matrix
# A = estimated mixing matrix
# S = estimated source matrix

library(data.table)
library(tidyverse)
library(fastICA)
library(fdrtool)
library(nortest)
library(factoextra)


setwd("/sc/arion/projects/psychgen/alanna/confound_comparisons/results/modules_ica/muscle_basic/")

# use number of ICs that reconstruct 95% of the variance as calculated by PCA
ge.mat <- fread("/sc/arion/projects/psychgen/alanna/confound_comparisons/data/expr/basic_normalized_expr/muscleskeletal_v8_basic_norm_expr.txt") %>% column_to_rownames("V1")
#ge.mat <- as.data.frame(t(ge.mat))
ge.mat <- as.data.frame(ge.mat)

PC <- prcomp(t(ge.mat), center=TRUE, scale.=FALSE)

project.pca.proportionvariances <- ((PC$sdev^2) / (sum(PC$sdev^2)))*100

n_PC = as.data.frame(table(cumsum(project.pca.proportionvariances)>95))[1,]$Freq+1
# n_PCs = n_FALSE + 1

# center and standardize expression of each gene across individuals before ICA
ge.mat <- scale(t(ge.mat),center=TRUE,scale=TRUE)
ge.mat <- as.data.frame(t(ge.mat))

#n = as.numeric(c("50","100","150","200","250","300","350","400"))

n=as.numeric(n_PC)
cuts = as.numeric(c("10e-2","10e-3","10e-4","10e-5","10e-6","10e-7","10e-8","10e-9","10e-10"))


for (i in 1:length(n)) {
	ica <- fastICA(ge.mat, n[i], fun="logcosh")
	for (k in 1:length(cuts)) {
		for (j in seq(1,n[i],1)) {
		dat <- ica$S[,j]
		qvals <- fdrtool(dat, cutoff.method="fndr")$qval
		sig <- (qvals < cuts[k])
		write.table(names(dat[sig]), paste0("genes_modnum",j,"_totalmods_",n[i],"_cutoff_",cuts[k]), quote=FALSE, row.names=FALSE, col.names=FALSE)
		} 
	}
write.table(ica$A, paste0(n[i], "_modules_ica_mixing_matrix.txt"), row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(ica$S, paste0(n[i], "_modules_ica_source_matrix.txt"), row.names=TRUE, col.names=TRUE, quote=FALSE)
}










