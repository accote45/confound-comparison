
# script for heart tissue dataset provided as example


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
library(dplyr)
library(sva)
library(pamr)




### normalize read counts btw samples using TMM in edgeR

expr.count <- fread("heart_left_ventricle_only_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
expr.count <- as.data.frame(expr.count) %>% column_to_rownames("V1")

expr.tpm <- fread("heart_left_ventricle_only_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
expr.tpm <- as.data.frame(expr.tpm) %>% column_to_rownames("V1")

index <- colnames(expr.count)

demos.attribs <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header=TRUE, fill=TRUE, sep="\t",quote = "",stringsAsFactors=FALSE)
demos.attribs <- demos.attribs[demos.attribs$SAMPID %in% index,]
demos.attribs$SampleID <- demos.attribs$SAMPID
demos.attribs <- separate(demos.attribs, SAMPID, c("V1","V2","V3","V4","V5"))
demos.attribs$SUBJID <- paste0(demos.attribs$V1,"-",demos.attribs$V2)

demos.phenos <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", fill=TRUE, header=TRUE)

covs <- read.table("Heart_Left_Ventricle.v8.covariates.txt", header=TRUE) %>% column_to_rownames("ID")
covs <- as.data.frame(t(covs)) %>% rownames_to_column("ID")
covs$SUBJID <- gsub("\\.","\\-",covs$ID)

temp <- merge(demos.attribs,demos.phenos,by="SUBJID")
demos.master <- merge(temp,covs,by="SUBJID")

#### find missing values
colSums(is.na(demos.master))


expr_fin <- expr.count[,colnames(expr.count) %in% index]

d <- DGEList(counts=expr_fin)
d <- calcNormFactors(d, method="TMM")

### filter genes based on threshold of >0.1 TPM in >/= 20% of samples and >/= 6 reads (unnormalized) in >/= 20% of samples

filt.genes.count <- rownames(expr.count[rowSums(expr.count>=6) >= 0.2*ncol(expr.count),])
filt.genes.tpm <- rownames(expr.tpm[rowSums(expr.tpm>=0.1) >= 0.2*ncol(expr.tpm),])
isexpr <- intersect(filt.genes.tpm, filt.genes.count)

d <- d[isexpr,,keep.lib.sizes=FALSE]


### Winsorize counts (Gene outliers in samples)
# Set gene counts in specific samples that are deviating 3 sd from other samples to 3 SD limit
wins.mat = apply(d$counts, 1, function(x){
  mn = mean(x, na.rm = T)
  std.dev = sd(x, na.rm = T)
  x[x < (mn-3*std.dev)] <- (mn-3*std.dev)
  x[x > (mn+3*std.dev)] <- (mn+3*std.dev)
  return(x)
}) %>% t

d$counts=wins.mat


## correlation between covariates 

# subset for only variables without na values
cor.table <- demos.master[,!(colnames(demos.master) %in% c("SMNUMGPS", "SM550NRM", "SM350NRM", "SMMNCPB", "SMMNCV", "SMCGLGTH", "SMGAPPCT", "SMNUM5CD", "DTHHRDY","sex","SUBJID","ID","SampleID"))]
cor.table <- cor.table[,!(colnames(cor.table) %like% "Inferred")]

# identify variables with no variance, remove
cor.table.fin <- Filter(function(x)(length(unique(x))>1), cor.table)
colnames(cor.table.fin)


form <- ~ SMATSSCR+	SMCENTER+	SMPTHNTS+	SMRIN+	SMTSISCH+	SMTSPAX+	SMNABTCH+	SMNABTCHT+	SMNABTCHD+	SMGEBTCH+	SMGEBTCHD+	SME2MPRT+	SMCHMPRS+	SMNTRART+	SMMAPRT+	SMEXNCRT+	SMGNSDTC+	SME1MMRT+	SMSFLGTH+	SMMPPD+	SMNTERRT+	SMRRNANM+	SMRDTTL+	SMVQCFL+	SMTRSCPT+	SMMPPDPR+	SMNTRNRT+	SMMPUNRT+	SMEXPEFF+	SMMPPDUN+	SME2MMRT+	SME2ANTI+	SMALTALG+	SME2SNSE+	SMMFLGTH+	SME1ANTI+	SMSPLTRD+	SMBSMMRT+	SME1SNSE+	SME1PCTS+	SMRRNART+	SME1MPRT+	SME2PCTS+	SEX+	AGE+	pcr+	platform
# compute canonical correlation analysis between all pairs of variables, output = rho/ sum(rho) (range of values = 0 to 1)
C = canCorPairs(form, demos.master)

pdf("allcovar_cca_gtex_heart.pdf")
plotCorrMatrix(C)
dev.off()

# remove collinear variables, repeat

C <- as.data.frame(C)
write.table(C,"/hpc/users/cotea02/heart_corrs.txt", quote=FALSE)


### variancePartition before covariate adjustment of expression matrix
demos.master$SEX <- as.character(demos.master$SEX)
demos.fin <- demos.master[demos.master$SampleID %in% index,]



# fit linear mixed model on gene expression (regression model fit to single gene), returns variance explained by each variable after correcting for other variables


# specify variables to be included in voom() estimates of uncertainty
design <- model.matrix(~ SMRIN + SMTSISCH, demos.fin)
vobjGenes <- voom(d, design)

form <- ~ SMRIN+	SMTSISCH+	SMTSPAX+	SMCHMPRS+	SMNTRART+	SMEXNCRT+	SME1MMRT+	SMVQCFL+	SMTRSCPT+	SMALTALG+	SMMFLGTH+	SMSPLTRD+	SMBSMMRT+	SME1SNSE+	SME1PCTS+	SMRRNART+	SME1MPRT+	SME2PCTS
cl <- makeCluster(8)
registerDoParallel(cl)

# in case of heart, had to use coefficient threshold of >0.88 --> so removed SMNTRNRT as well, otherwise variancePartition gave error of strong collinearity 


varPart <- fitExtractVarPartModel(vobjGenes, form, demos.fin)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )


# Bar plot of variance fractions for the first 10 genes
pdf("heart_gtex_barplot_variancefractions_first10genes_unadjustedGE.pdf")
plotPercentBars( vp[1:10,] )
dev.off()

# violin plot of contribution of each variable to total variance
pdf("heart_gtex_violin_varcontribution_unadjustedGE.pdf")
plotVarPart( vp )
dev.off()


write.table(varPart, "heart_varpart_output.txt", quote=FALSE)



# keep covariates that explain >/=1% of expr variation in >/=10% of genes

colSums(varPart>=0.01) >= 0.1*nrow(varPart)

colSums(dat>=0.01) >= 0.1*nrow(dat)


### Normalisation
 
# design matrix
design.adj <- model.matrix(~0 + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + SMRIN + SMTSISCH + SMNTRART + SMEXNCRT + SMVQCFL + SMTRSCPT + SMALTALG + SMMFLGTH + SMSPLTRD + SMBSMMRT + SME1SNSE + SME1PCTS + SMRRNART + SME1MPRT + SME2PCTS , demos.fin)


# Estimate voom weights
#tmp = y$counts
#tmp[is.na(tmp)] = 0
VOOM.GENE_EXPRESSION = voom(d, design=design.adj, plot=F)

# Fit linear model using new weights and new design
ADJUSTED.FIT = lmFit(VOOM.GENE_EXPRESSION)
  
# Residuals after normalisation
RESIDUAL.GENE_EXPRESSION = residuals.MArrayLM(ADJUSTED.FIT, VOOM.GENE_EXPRESSION)
write.table(RESIDUAL.GENE_EXPRESSION, "heart_knowncovar_wAnc_adj.txt", quote=FALSE)
