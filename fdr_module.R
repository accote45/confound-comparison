


###########################################
### WGCNA analysis

library(tidyverse)
library(plyr)
library(biomaRt)
library(purrr)
library(gprofiler2)
library(data.table)
library(biomaRt)
library(httr)
library(jsonlite)
library(GeneOverlap)
library(reshape2)
library(pracma)
library(qvalue)

coexpr=args[1]

setwd(paste0("/modules_wgcna/",coexpr))

listfiles <- list.files(pattern=glob2rx("unsigned_pearson_min_mod_size10module*"))
listfiles <- listfiles[!(listfiles %in% c("unsigned_pearson_min_mod_size10module_assigned_master.txt","unsigned_pearson_min_mod_size10module_eigengenes_10min.txt","unsigned_pearson_min_mod_size10module_grey"))]

# read in coexpression modules
coexpr_mods <- list()
for (k in 1:length(listfiles)){
	coexpr_mods[[k]] <- read.table(text=readLines(listfiles[k]), header=FALSE, quote="")
	coexpr_mods[[k]] <- gsub("\\..*","",coexpr_mods[[k]]$V1) 	
}
names(coexpr_mods) <- listfiles

dat <- fread("expression.txt",stringsAsFactors=FALSE) %>% column_to_rownames(var = "V1")
#dat <- t(dat)
rownames(dat) <- gsub("\\..*","",rownames(dat))

ref <- read.table("/somekh_evaluation/tissue_network_evaluation/intestine_gold_standard2.txt")

library(biomaRt)
# Define biomart object
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# Query biomart
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "gene_biotype"),
                       filters = "ensembl_gene_id", values = rownames(dat),
                       mart = mart, uniqueRows = TRUE)

new <- ref 
new[] <- lapply(ref[,1:2], function(x) Ensemble2HGNC$ensembl_gene_id[match(x, Ensemble2HGNC$entrezgene_id)])
new$V3 <- ref$V3


# subset gold.standard for genes present in input data files
index <- rownames(dat)
new.fin <- new[(new$V1 %in% index & new$V2 %in% index),]

# count/document number of true pos and true neg pairs used for evaluation of this expr dataset
dim(new.fin[new.fin$V3=="1",])
dim(new.fin[new.fin$V3=="0",])

# subset new for true pos cases 
# if new$V1 and new$V2 in module, print ID in new vector (tp)

tp = c()
new.pos <- new[new$V3=="1",]

for(j in 1:length(coexpr_mods))
{
	count <- ifelse(new.pos$V1 %in% coexpr_mods[[j]] & new.pos$V2 %in% coexpr_mods[[j]],"1","0")
	tp <- c(tp,count)
}
tp.count <- as.data.frame(table(tp))[2,]$Freq


# subset new for true negative cases 
# if new$V1 and new$V2 in module, print ID in new vector (tn)

fp = c()
new.neg <- new[new$V3=="0",]

for(j in 1:length(coexpr_mods))
{	count <- ifelse(new.neg$V1 %in% coexpr_mods[[j]] & new.neg$V2 %in% coexpr_mods[[j]],"1","0") # in network and in true negative class
	fp <- c(fp,count)
}

fp.count <- as.data.frame(table(fp))[2,]$Freq
tn.count <- 5000-as.data.frame(table(fp))[2,]$Freq
fn.count <- 5000-tp.count

# calculate FP, TP, FP
precision <- tp.count/(tp.count+fp.count)
recall <- tp.count/(tp.count+fn.count)
fdr <- fp.count/(tp.count+ fp.count)
fnr <- fn.count/(fn.count+tp.count)

cat(c(fdr,fnr,precision,recall,tp.count,fp.count,tn.count,fn.count))




