


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



expr.count <- fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",skip=2, header=TRUE)
expr.count <- as.data.frame(expr.count) %>% column_to_rownames("Name")

expr.tpm <- fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", skip=2, header=TRUE, sep="\t")
expr.tpm <- as.data.frame(expr.tpm) %>% column_to_rownames("Name")

demos.attribs <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header=TRUE, fill=TRUE, sep="\t",quote = "",stringsAsFactors=FALSE)
#demos.phenos <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")


### subset expr data for tissue of interest (only complete once, save files)

# create expr sample ID
demos.attribs$SampleID <- demos.attribs$SAMPID
demos.attribs <- separate(demos.attribs, SAMPID, c("V1","V2","V3","V4","V5"))
demos.attribs$EXPRID <- paste0(demos.attribs$V1,"-",demos.attribs$V2)

# subset expression dataset for samples with WGS + RNA seq for tissue of interest

index <- demos.attribs[demos.attribs$SMAFRZE=="WGS",]$EXPRID
demos.attribs.withgen <- demos.attribs[demos.attribs$EXPRID %in% index,]

index.expr <- demos.attribs.withgen[demos.attribs.withgen$SMTSD %in% c("Spleen") & demos.attribs.withgen$SMAFRZE=="RNASEQ",]$SampleID

expr.count.tissue <- expr.count[,colnames(expr.count) %in% index.expr]
write.table(expr.count.tissue, "spleen_only_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", quote=FALSE)


expr.tpm.tissue <- expr.tpm[,colnames(expr.tpm) %in% index.expr]
write.table(expr.tpm.tissue, "spleen_only_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", quote=FALSE)






####################################################################################################################################


### normalize read counts btw samples using TMM in edgeR

expr.count <- fread("spleen_only_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
expr.count <- as.data.frame(expr.count) %>% column_to_rownames("V1")

expr.tpm <- fread("spleen_only_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
expr.tpm <- as.data.frame(expr.tpm) %>% column_to_rownames("V1")

demos.attribs <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header=TRUE, fill=TRUE, sep="\t",quote = "",stringsAsFactors=FALSE)

d <- DGEList(counts=expr.count)
d <- calcNormFactors(d, method="TMM")

# return TMM normalized log2CPM
tmm.cpm <- cpm(d, normalized.lib.sizes=TRUE, log=TRUE)




### filter genes based on threshold of >0.1 TPM in >/= 20% of samples and >/= 6 reads (unnormalized) in >/= 20% of samples

filt.genes.count <- rownames(expr.count[rowSums(expr.count>=6) >= 0.2*ncol(expr.count),])
filt.genes.tpm <- rownames(expr.tpm[rowSums(expr.tpm>=0.1) >= 0.2*ncol(expr.tpm),])
isexpr <- intersect(filt.genes.tpm, filt.genes.count)

tmm.cpm.expr <- tmm.cpm[rownames(tmm.cpm) %in% isexpr,]




### Winsorize counts (Gene outliers in samples)
# Set gene counts in specific samples that are deviating 3 sd from other samples to 3 SD limit
cpm.wins.mat = apply(tmm.cpm.expr, 1, function(x){
  mn = mean(x, na.rm = T)
  std.dev = sd(x, na.rm = T)
  x[x < (mn-3*std.dev)] <- (mn-3*std.dev)
  x[x > (mn+3*std.dev)] <- (mn+3*std.dev)
  return(x)
}) %>% t


write.table(cpm.wins.mat,"spleen_v8_basic_norm_expr.txt",quote=FALSE, col.names=TRUE, row.names=TRUE)

