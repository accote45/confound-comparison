

library(peer)
library(tidyverse)

expr <- args[1]
tissue <- args[2]

#build model
model=PEER()

# set max number of factors
PEER_setNk(model,30)

# set expression data
PEER_setPhenoMean(model,as.matrix(expr))

# add factor to account for mean expression
#PEER_setAdd_mean(model, TRUE)

# train the model (100 iterations should be sufficient to reach convergence on most datasets)
PEER_update(model)

#plot posterior variance of the factor weights and convergence diagnostics - if there is an elbow, consider only including the more relevant factors
#pdf('test.pdf')
#PEER_plotModel(model)
#dev.off()

# get factors (posterior mean of the inferred confounders)
factors <- as.data.frame(PEER_getX(model))
rownames(factors) <- rownames(expr)
write.table(factors, paste0(tissue,"_peer_factors.txt"),quote=FALSE, col.names=TRUE, row.names=TRUE)

# get weights of inferred confounders
weights = PEER_getW(model)
rownames(weights)<- colnames(expr)
write.table(weights, paste0(tissue,"_peer_weights.txt"),quote=FALSE, col.names=TRUE, row.names=TRUE)

# precision (inverse variance) of the weights
precision = PEER_getAlpha(model)
write.table(precision, paste0(tissue,"_peer_precision.txt"),quote=FALSE)

# get residuals
residuals <- PEER_getResiduals(model)
rownames(residuals) <- rownames(expr)
colnames(residuals) <- colnames(expr)
write.table(residuals, paste0(tissue,"_peer_adj.txt"),quote=FALSE, col.names=TRUE, row.names=TRUE)
