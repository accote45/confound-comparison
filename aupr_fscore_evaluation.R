


#############################
# calculate AUPR

library(tidyverse)
library(dplyr)
library(data.table)
library(pROC)
library(matrixStats)

dat <- args[1]
#dat <- t(dat)
rownames(dat) <- gsub("\\..*","",rownames(dat))

ref <- read.table("/tissue_network_evaluation/blood_gold_standard2.txt")

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
#new$V1 <- gsub("\\..*","",new$V1)
#new$V2 <- gsub("\\..*","",new$V2)


# subset gold.standard for genes present in input data files
index <- rownames(dat)
new.fin <- new[(new$V1 %in% index & new$V2 %in% index),]

# count/document number of true pos and true neg pairs used for evaluation of this expr dataset
dim(new.fin[new.fin$V3=="1",])
dim(new.fin[new.fin$V3=="0",])


source("PRROC.R")


# calculate correlations for true pos/neg gene pairs
  bins.cors.df                  = new.fin
  #bins.cors.df$binAll.raw.r     = NA#add a zero column
  bins.cors.df$binAll.raw.pval  = NA

  for(j in 1:nrow(new.fin))
  {
    gene1.ensmbl = new.fin$V1[j]
    gene2.ensmbl = new.fin$V2[j]
    
    if( (gene1.ensmbl %in% rownames(dat)) & (gene2.ensmbl %in% rownames(dat)))
    {
      gene1.bin.all.raw.vec   = gene2.bin.all.raw.vec = 0
      gene1.bin.all.raw.vec   = dat[gene1.ensmbl,]
      gene2.bin.all.raw.vec   = dat[gene2.ensmbl,]
      result.all.raw          = stats::cor.test(as.numeric(gene1.bin.all.raw.vec), as.numeric(gene2.bin.all.raw.vec), method = "pearson")
      
      #bins.cors.df[j,"binAll.raw.r"]    = result.all.raw$estimate
      bins.cors.df[j,"binAll.raw.pval"] = result.all.raw$p.value
      #bins.cors.df[j,"binAll.r"]        = result.all$estimate 
    }else{
      print(paste0(j," is index of not found gene in input.edata"))
      j = j+1
    }
  }


  bins.cors.df$BH <- p.adjust(bins.cors.df[,4], method = "BH")#multiple correction for each column
  bins.cors.df$log <- -log(bins.cors.df$BH,10)


 temp <- bins.cors.df$log
 temp[is.infinite(temp)]   = 0
 max.val                 = max(temp)
 bins.cors.df$log[is.infinite(bins.cors.df$log)]<-(max.val+10)#change Inf to max value

fg <- bins.cors.df[bins.cors.df$V3=="1",]$log # pvals for each true positive gene pair
bg <- bins.cors.df[bins.cors.df$V3=="0",]$log # pvals for each true negaitve gene pair

pr <- pr.curve(scores.class0=fg, scores.class1 = bg, curve=T,rand.compute=TRUE)

pr$rand
pr$auc.integral

pdf("aupr.pdf")
plot(pr)
dev.off()













###################
### Calculated F-score

# harmonic mean of precision and recall at a given cutoff
# will test F-score at different thresholds for gene-gene pair (pearson coefficient 0.3-0.7)

library(tidyverse)
library(dplyr)
library(data.table)

dat <- args[1]
rownames(dat) <- gsub("\\..*","",rownames(dat))

ref <- read.table("/sc/arion/projects/psychgen/alanna/confound_comparisons/results/somekh_evaluation/tissue_network_evaluation/adipose_gold_standard.txt")

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


# calculate correlations, keep with absolute pearson > 0.5

dat <- as.matrix(dat)

  bins.cors.df                  = new
  #bins.cors.df$binAll.raw.r     = NA#add a zero column
  bins.cors.df$binAll.corr  = NA

for(j in 1:nrow(new))
  {
    gene1.ensmbl = new$V1[j]
    gene2.ensmbl = new$V2[j]
    
    if( (gene1.ensmbl %in% rownames(dat)) & (gene2.ensmbl %in% rownames(dat)))
    {
      gene1.bin.all.raw.vec   = gene2.bin.all.raw.vec = 0
      gene1.bin.all.raw.vec   = dat[gene1.ensmbl,]
      gene2.bin.all.raw.vec   = dat[gene2.ensmbl,]
      result.all.raw          = stats::cor.test(gene1.bin.all.raw.vec, gene2.bin.all.raw.vec, method = "pearson")
      
      #bins.cors.df[j,"binAll.raw.r"]    = result.all.raw$estimate
      bins.cors.df[j,"binAll.corr"] = result.all.raw$estimate
  }
}

# calculate FP, TP, FP

total.pos=5000
total.neg=5000

master <- data.frame(matrix(NA,nrow=9,ncol=8))
colnames(master) <- c("tp","fn","fp","precision","recall","fdr","fnr","f_score")

cutoff <- c(0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7)
rownames(master) <- cutoff
for(i in 1:length(cutoff))
{
  master$tp[i] <- nrow(bins.cors.df[(bins.cors.df$V3=="1") & (abs(bins.cors.df$binAll.corr)>cutoff[i]),])
  master$fn[i] <- nrow(bins.cors.df[(bins.cors.df$V3=="1") & (abs(bins.cors.df$binAll.corr)<cutoff[i]),]) # not in network, in true positive class
  master$fp[i] <- nrow(bins.cors.df[(bins.cors.df$V3=="0") & (abs(bins.cors.df$binAll.corr)>cutoff[i]),]) # in network, not in true positive class
  master$precision[i] <- master$tp[i]/(master$tp[i]+master$fp[i])
  master$recall[i] <- master$tp[i]/(master$tp[i]+master$fn[i])
  master$fdr[i] <- master$fp[i]/(master$tp[i]+ master$fp[i])
  master$fnr[i] <- master$fn[i]/(master$fn[i]+master$tp[i])
  master$f_score[i] <- 2*(master$precision[i]*master$recall[i])/(master$precision[i]+master$recall[i])
}

write.table(master,'f_score.txt',quote=F)





