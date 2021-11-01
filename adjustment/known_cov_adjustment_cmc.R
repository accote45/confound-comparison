

####################################################################################
######### Control RNA seq for covariates, create matrix of residualized GE ##########
#####################################################################################

# Adjust raw CMC DLPFC RNA seq for covariates through residualizing
	# Voom normalization and gene filtering
	# Define covs through correlation of covs with PCA of voom-normalized expression
	# Calculate residuals to remove confounds effect from expression matrix

setwd("/data")

library(tidyverse)
library(biomaRt)
library(plyr)
library(limma)
library(Glimma)
library(edgeR)
library('variancePartition')
library('doParallel')
library(data.table)
library(corrplot)
library(RColorBrewer)

library(MatrixEQTL)

COUNT <- fread("expr_files/txicounts_+HBCC.csv.gz",header=T,sep=",",stringsAsFactors=F)
demos <- read.csv("demos/demos.master.EA.CMC1.amanda.ancestry.only.csv",header=T,stringsAsFactors=F) 

# Create RIN^2 and clustered lib batch covariates
demos$RIN2 <- demos$rnaSeq_isolation.RIN^2

cBatchToBatch=as.matrix(list("0"=c("11/11/13","11/26/13"), "A"=c(17,18,10,"8/28/13"), "B"=c(14,7,11,6,16,25), "C"=c(5,1,3), "D"=c(2,28,8,12,4), "E"=c(20,24,26,27,9), "F"=c(13,15,22,21,23), "G"=c("10/15/13",19), "H"=c("10/9/13")))
key <- as.data.frame(cBatchToBatch) %>% rownames_to_column("clusterLIB")
key$Library_Batch.2 <- sapply(key$V1, paste, collapse=",")
key <- key %>% separate_rows(Library_Batch.2, sep = "," )
key$V1 <- NULL
demos$Library_Batch.2 <- substring(demos$Library_Batch, 3)

demos <- merge(demos, key, by.x="Library_Batch.2")


# Define biomart object
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",mirror="asia")


# Query biomart
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "percentage_gene_gc_content", "gene_biotype", "chromosome_name", "start_position", "end_position"),
                       filters = "ensembl_gene_id", values = COUNT$V1,
                       mart = mart, uniqueRows = TRUE)

#Ensemble2HGNC <- read.table("/sc/hydra/projects/psychgen/alanna/biomart_export_GRCh38.p13.txt", sep=",", header=TRUE)
Ensemble2HGNC.dedup <- Ensemble2HGNC[!duplicated(Ensemble2HGNC$ensembl_gene_id),]


COUNT$ensembl_gene_id <- COUNT$V1
COUNT = join(COUNT, Ensemble2HGNC.dedup, by="ensembl_gene_id", type="left") %>% filter(chromosome_name != "X" & chromosome_name != "Y" & chromosome_name !="MT")
COUNT <- as.data.frame(COUNT)

# Keep only individuals with info in demos file
index <- c(demos$RNAseq.Sample_RNA_ID, "ensembl_gene_id")
count.trim <- COUNT[,(names(COUNT) %in% index)] %>% column_to_rownames("ensembl_gene_id")

## To remove participants missing one or more covars
covars <- c("Genotypes.Genotyping_Sample_ID","RNAseq.Sample_RNA_ID",'Dx', 'Sex', 'Age_of_Death', 'Institution', 'PMI_.in_hours.', "rnaSeq_isolation.RIN", "RIN2", "EV.1", "EV.2", "EV.3", "EV.4", "EV.5", "clusterLIB")
demos.fin <- demos[,covars] %>% na.omit()
keep <- demos.fin$RNAseq.Sample_RNA_ID
rownames(demos.fin) <- NULL
count.master <- count.trim[ , (names(count.trim) %in% keep)]

### Explore correlation among covariates
covar_corr <- cor(demos.fin[,unlist(lapply(demos.fin, is.numeric))], method='s', use='pairwise.complete.obs')     
pdf("/figures/contcovar_corr_cmc1_euro_dlpfc.pdf")
corrplot(as.matrix(covar_corr), method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 4, rect.col = "black", rect.lwd = 5,cl.pos = "b", tl.col = "indianred4", tl.cex = 1.0, cl.cex = 1.0, addCoef.col = "white", number.digits = 2, number.cex = 0.50, col = colorRampPalette(c("darkred","white","midnightblue"))(100))
dev.off()

form <- ~ Dx + Sex + Age_of_Death + Institution + PMI_.in_hours. + rnaSeq_isolation.RIN + RIN2 + EV.1 + EV.2 + EV.3 + EV.4 + EV.5 + clusterLIB

# compute canonical correlation analysis between all pairs of variables, output = rho/ sum(rho) (range of values = 0 to 1)
C = canCorPairs(form, demos.fin)

pdf("/figures/allcovar_cca_cmc1_euro_dlpfc.pdf")
plotCorrMatrix(C)
dev.off()


### Normalisation

# Sort files, make DGEList object
demos.sort <- demos.fin[order(demos.fin$RNAseq.Sample_RNA_ID),] 
count.sort <- count.master[ , order(names(count.master))]

y = DGEList( count.sort, genes=rownames(count.sort) )
isexpr = rowSums(cpm(count.sort)>1) >= 0.5*ncol(count.sort)
y <- y[isexpr,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
dim(y)

sizeFactors.subset <- as.data.frame(t(y$samples))[2,]
index <- rownames(y$counts)
uCovar <- Ensemble2HGNC.dedup[(Ensemble2HGNC.dedup$Gene.stable.ID %in% index),]
rownames(uCovar) <- uCovar$Gene.stable.ID
uCovar$length <- uCovar$Gene.end..bp. - uCovar$Gene.start..bp.
uCovar <- uCovar[,c(5,9)]
uCovar<-uCovar[sort(rownames(uCovar)),]


# CQN
library(cqn)
### Library Normalisation
# Compute offset for gene length and gc content
cqn.expr = cqn(y$counts, lengths=uCovar$length, x=uCovar$Gene...GC.content, lengthMethod = "smooth", verbose=TRUE)

# normalized expr matrix (values on log2RPKM scale)
# RPKM = (counts / (total reads / 1,000,000)) / length of gene in KB
cqn.norm.expr <- cqn.expr$y + cqn.expr$offset


y = DGEList( y$counts, genes=uCovar)
y$offset <- cqn.expr$glm.offset
dim(y)



### Removal of outlier samples 
# Find principal components of expression to plot
PC <- prcomp(t(cqn.norm.expr), scale.=T, center = T)
# Plot first 2 PCs
plotdata <- data.frame(RNAseq.Sample_RNA_ID=rownames(PC$x), 
                       PC1=PC$x[,1], 
                       PC2=PC$x[,2])
plotdata <- merge(plotdata, demos.sort, by="RNAseq.Sample_RNA_ID", all.y=TRUE)

pdf("/figures/PCAplot_cmc1_euro_dlpfc_unadjusted.pdf")
p <- ggplot(plotdata, aes(x=PC1, y=PC2))
p <- p + geom_point(aes(colour=factor(Institution), fill=factor(Institution), shape=Dx, size=Age_of_Death)) 
p <- p + scale_shape_manual(values=c(3,22,21,24))
p <- p + theme_bw() %+replace% theme(legend.position="right")
p
dev.off()

# For CMC1 cohort only, removed single outlier based on inspection of pc plot
plotdata <- plotdata[order(plotdata$PC2),]
tail(plotdata)

indremove <- c("MSSM_RNA_PFC_139")
cqn.norm.expr = cqn.norm.expr[,!colnames(cqn.norm.expr) %in% indremove]
y = y[,!colnames(y) %in% indremove]
demos.sort = demos.sort[!demos.sort$RNAseq.Sample_RNA_ID %in% indremove,]

# check for outlier samples by inter array correlation (samples with IAC < 3SDs below mean IAC for the dataset will be removed) ### source = https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/HumanBrainTranscriptome/Identification%20and%20Removal%20of%20Outlier%20Samples%20-%20Illumina.pdf
IAC=cor(cqn.norm.expr,use="p") 
library(WGCNA)

png("/figures/iac_histogram_cmc1.png")
hist(IAC,sub=paste("mean=",format(mean(IAC[upper.tri(IAC)]),digits=3))) 
dev.off()

sampleTree = hclust(as.dist(1-IAC), method = "average");
png("/figures/iacbased_outlier_check_cmc1.png")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
abline(h=0.1,col="red") 
dev.off()

groups.2 <- cutree(sampleTree,2)
table(groups.2)
groups.2[order(groups.2)]

meanIAC=apply(IAC,2,mean)
sdCorr=sd(meanIAC)
numbersd=(meanIAC-mean(meanIAC))/sdCorr
png('/figures/iac_sd_graph_cmc1.png')
plot(numbersd)
abline(h=-3)
dev.off()

sdout=-3
outliers=dimnames(cqn.norm.expr)[[2]][numbersd<sdout]
outliers

removevec=c("MSSM_RNA_PFC_139", "MSSM_RNA_PFC_133", "MSSM_RNA_PFC_163", "MSSM_RNA_PFC_299", "MSSM_RNA_PFC_304", "MSSM_RNA_PFC_60", "PENN_RNA_PFC_56", "PENN_RNA_PFC_71")
overlap1=is.element(dimnames(cqn.norm.expr)[[2]],removevec)
datrest2=cqn.norm.expr[,!overlap1]
dim(datrest2)
# [1] 20510 484
IAC=cor(datrest2,use="p")
png("/figures/iac_histogram_afteroutlierremoval_cmc1.png")
hist(IAC,sub=paste("mean=",format(mean(IAC[upper.tri(IAC)]),digits=3))) 
dev.off()

cluster1=hclust(as.dist(1-IAC),method="average")
png("/figures/iacbased_outlier_check_round2_cmc1.png")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(cluster1, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
abline(h=0.1,col="red") 
dev.off()

indremove=c("MSSM_RNA_PFC_139", "MSSM_RNA_PFC_133", "MSSM_RNA_PFC_163", "MSSM_RNA_PFC_299", "MSSM_RNA_PFC_304", "MSSM_RNA_PFC_60", "PENN_RNA_PFC_56", "PENN_RNA_PFC_71")
cqn.norm.expr = cqn.norm.expr[,!colnames(cqn.norm.expr) %in% indremove]
y = y[,!colnames(y) %in% indremove]
demos.sort = demos.sort[!demos.sort$RNAseq.Sample_RNA_ID %in% indremove,]

### Winsorize counts (Gene outliers in samples)
# Set gene counts in specific samples that are deviating 3 sd from other samples to 3 SD limit
raw.wins.mat = apply(cqn.norm.expr, 1, function(x){
  mn = mean(x, na.rm = T)
  std.dev = sd(x, na.rm = T)
  x[x < (mn-3*std.dev)] <- (mn-3*std.dev)
  x[x > (mn+3*std.dev)] <- (mn+3*std.dev)
  return(x)
}) %>% t

cqn.norm.expr = raw.wins.mat

raw.wins.mat = apply(y$counts, 1, function(x){
  mn = mean(x, na.rm = T)
  std.dev = sd(x, na.rm = T)
  x[x < (mn-3*std.dev)] <- (mn-3*std.dev)
  x[x > (mn+3*std.dev)] <- (mn+3*std.dev)
  return(x)
}) %>% t

y$counts = raw.wins.mat



### variancePartition before covariate adjustment of expression matrix

cl <- makeCluster(4)
registerDoParallel(cl)

# fit linear mixed model on gene expression (regression model fit to single gene), returns variance explained by each variable after correcting for other variables
form <- ~ (1|Dx) + (1|Sex) + Age_of_Death + (1|Institution) + PMI_.in_hours. + rnaSeq_isolation.RIN + RIN2 + EV.1 + EV.2 + EV.3 + EV.4 + EV.5 + (1|clusterLIB)
varPart <- fitExtractVarPartModel(cqn.norm.expr, form, demos.sort)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )

# Bar plot of variance fractions for the first 10 genes
pdf("/figures/cmc1_euro_barplot_variancefractions_first10genes_unadjustedGE.pdf")
plotPercentBars( vp[1:10,] )
dev.off()

# violin plot of contribution of each variable to total variance
pdf("/figures/cmc1_euro_violin_varcontribution_unadjustedGE.pdf")
plotVarPart( vp )
dev.off()

write.table(cqn.norm.expr, "expr_files/quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt")
write.table(cqn.norm.expr, "/quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt")


### Normalisation

# Estimate voom weights
#tmp = y$counts
#tmp[is.na(tmp)] = 0
VOOM.GENE_EXPRESSION = voom(y, plot=F)


### (FOR CHECK) PCA of adjusted GE matrix 
PC <- prcomp(VOOM.GENE_EXPRESSION, scale.=T, center = T)
# Plot first 2 PCs
plotdata <- data.frame(RNAseq.Sample_RNA_ID=rownames(PC$rotation), 
                       PC1=PC$rotation[,1],
                       PC2=PC$rotation[,2])
plotdata <- merge(plotdata, demos.sort, by="RNAseq.Sample_RNA_ID", all.y=TRUE)

pdf("/figures/PCAplot_cmc1_euro_dlpfc_no_covaradj.pdf")
p <- ggplot(plotdata, aes(x=PC1, y=PC2))
#p <- p + geom_point(aes(colour=factor(Institution), fill=factor(Institution), shape=Dx, size=Age_of_Death)) 
p <- p + geom_point(aes(colour=factor(Institution), fill=factor(Institution), size=Age_of_Death)) 
p <- p + scale_shape_manual(values=c(3,22,21,24)) + scale_colour_manual(values=c("black", "black", "black", "black"))
p <- p + theme_bw() %+replace% theme(legend.position="right")
p
dev.off()

dat <- as.data.frame(VOOM.GENE_EXPRESSION)
write.table(dat, "expr_files/quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt")






### Normalisation
  
# design matrix
design.adj <- model.matrix(~0 + Dx + Sex + Age_of_Death + Institution + PMI_.in_hours. + rnaSeq_isolation.RIN + RIN2 + clusterLIB + EV.1 + EV.2 + EV.3 + EV.4 + EV.5, demos.sort)

# Estimate voom weights
#tmp = y$counts
#tmp[is.na(tmp)] = 0
VOOM.GENE_EXPRESSION = voom(y, design=design.adj, plot=F)
VOOM.GENE_EXPRESSION$E = cqn.norm.expr

# Fit linear model using new weights and new design
ADJUSTED.FIT = lmFit(VOOM.GENE_EXPRESSION)
  
# Residuals after normalisation
RESIDUAL.GENE_EXPRESSION = residuals.MArrayLM(ADJUSTED.FIT, cqn.norm.expr)





# variancePartition after covariate adjustment of expression matrix
cl <- makeCluster(4)
registerDoParallel(cl)

# fit linear mixed model on gene expression (regression model fit to single gene), returns variance explained by each variable after correcting for other variables
form <- ~ (1|Dx) + (1|Sex) + Age_of_Death + (1|Institution) + PMI_.in_hours. + rnaSeq_isolation.RIN + RIN2 + EV.1 + EV.2 + EV.3 + EV.4 + EV.5 + (1|clusterLIB)

varPart <- fitExtractVarPartModel(RESIDUAL.GENE_EXPRESSION, form, demos.sort)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )

# Bar plot of variance fractions for the first 10 genes
pdf("cmc1_euro_barplot_variancefractions_first10genes_adjustedGE_wAnc.pdf")
plotPercentBars( vp[1:10,] )
dev.off()

# violin plot of contribution of each variable to total variance
pdf("/cmc1_euro_violin_varcontribution_adjustedGE_wAnc.pdf")
plotVarPart( vp )
dev.off()

write.table(RESIDUAL.GENE_EXPRESSION, "/knowncovar_adj_quantnorm_outlierrem_winsorized_expression_cmc1_euro.txt")





