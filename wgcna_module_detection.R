


library(WGCNA)
library(data.table)
library(tidyverse)
library(dplyr)
library(plyr)

options(stringsAsFactors = FALSE);
allowWGCNAThreads(nThreads=5)

mod = 10

dat <- args[1]


gene.names=colnames(dat)

# check for samples and genes w/ too many missing entries, filter if necessary
gsg = goodSamplesGenes(dat, verbose = 3);
gsg$allOK

########### Network construction and module detection ##############

# Choose a set of soft-thresholding powers
powers = c(c(1:11), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(dat, dataIsExpr = TRUE, powerVector = powers, ,blockSize=26000, verbose = 5, networkType = "unsigned")
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
png(paste0("unsigned_pearson_min_mod_size",mod,"_wgcna_softthreshold_scaleind_results.png"))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
dev.off()

# Mean connectivity as a function of the soft-thresholding power
png(paste0("unsigned_pearson_min_mod_size",mod,"wgcna_softthreshold_meanconnect_results.png"))
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# specify power parameter
softPower=2

# calculate adjacencies 
adjacency = adjacency(dat, power = softPower);
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
colnames(TOM) =rownames(TOM) =gene.names

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

## clustering based on topological overlap
minModuleSize = mod;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)

## to plot module assignment under the dendrogram
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
png(paste0("unsigned_pearson_min_mod_size",mod,"wgcna_hclust_dendrogram_10min.png"))
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
dev.off()

## merge modules with similar expression profiles
#calc eigengenes
MEList = moduleEigengenes(dat, colors=dynamicColors)
MEs = MEList$eigengenes
#calculat dissimiliarity of module eigengenes
MEdiss = 1-cor(MEs);
#cluster MEs
METree = hclust(as.dist(MEdiss), method='average');
#plot result
png(paste0("unsigned_pearson_min_mod_size",mod,"wgcna_hclust_dendrogram_10min_mergedmodules.png"))
plot(METree, main = 'Clustering of module eigengenes')
MEDissThres=0.25
abline(h=MEDissThres,col='red')
dev.off()

merge=mergeCloseModules(dat,dynamicColors,cutHeight=MEDissThres,verbose=3)
#merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
png(paste0("unsigned_pearson_min_mod_size",mod,"wgcna_hclust_dendrogram_10min_mergedmodules.png"))
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dat, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

# To extract modules
module_colors= unique(mergedColors)
for (color in module_colors){
    module=gene.names[which(mergedColors==color)]
    write.table(module, paste0("unsigned_pearson_min_mod_size",mod,"module_",color), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

library(reshape2)
master <- melt(data.frame(gene.names,mergedColors))
write.table(master,paste0("unsigned_pearson_min_mod_size",mod,"module_assigned_master.txt"))


#eigengenes
write.table(as.data.frame(MEs), paste0("unsigned_pearson_min_mod_size",mod,"module_eigengenes_10min.txt"),quote=FALSE)
write.table(as.data.table(table(mergedColors)),paste0("unsigned_pearson_min_mod_size",mod,"numbergenes_permodule.txt"),quote=FALSE)

hubgenes <- chooseTopHubInEachModule(dat, mergedColors, omitColors = "grey",power = 2,type = "unsigned")
write.table(as.data.frame(hubgenes),paste0("unsigned_pearson_min_mod_size",mod,"hubgenes_rosmap.txt"),quote=FALSE)

connectivity <- intramodularConnectivity(adjacency, mergedColors,
              scaleByMax = FALSE)

write.table(as.data.frame(connectivity),paste0("unsigned_pearson_min_mod_size",mod,"connectivity.txt"), quote=FALSE)

# hierarchical clustering of MEs
dissimME=(1-t(cor(MEs, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
png(paste0("unsigned_pearson_min_mod_size",mod,"clustering_eigengenes.png"))
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")
dev.off()


#pearson correlation of gene with module eigengene
KME <- signedKME(dat, MEs, outputColumnName="kME")
write.table(as.data.frame(KME),paste0("unsigned_pearson_min_mod_size",mod,"KME.txt"), quote=FALSE)

# avg correlation of gene w/module eigengene per module
module_colors= unique(mergedColors)
n = length(module_colors)
avgs <- c()
name <- c()
for (i in 1:n){
restrictGenes=mergedColors==module_colors[i]
avgs[i] <- mean(KME[restrictGenes, paste0("kME",module_colors[i],sep="")])
name[i] <- module_colors[i]
}
names(avgs) <- name
kme.avg <- cbind(read.table(text=names(avgs)), avgs)
write.table(as.data.frame(KME),paste0("unsigned_pearson_min_mod_size",mod,"kme.avg.txt"), quote=FALSE)

# avg intramodular connectivity per module
module_colors= unique(mergedColors)
n = length(module_colors)
avgs <- c()
name <- c()
for (i in 1:n){
restrictGenes=mergedColors==module_colors[i]
avgs[i] <- mean(connectivity[restrictGenes, "kWithin"])
name[i] <- module_colors[i]
}
names(avgs) <- name
connect.avg <- cbind(read.table(text=names(avgs)), avgs)
write.table(as.data.frame(KME),paste0("unsigned_pearson_min_mod_size",mod,"connectivity.avg.txt"), quote=FALSE)



