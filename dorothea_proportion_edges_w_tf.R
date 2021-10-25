


### percentage of network edges involving a TF


library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(plyr)
library(WGCNA)
library(matrixStats)
library(beeswarm)

setwd("/scripts/edge_list/cutoff_0.3")

files2 <- list.files(pattern = "*edge.list", full.names=FALSE)
#files2 <- files2[grepl("intestine",files2)]

listOfFiles2 <- lapply(files2, function(x) fread(x,drop=1)) 

for (i in 1:length(listOfFiles2)){
	listOfFiles2[[i]]$V1 <- gsub("\\..*","",listOfFiles2[[i]]$V1)
	listOfFiles2[[i]]$V2 <- gsub("\\..*","",listOfFiles2[[i]]$V2)
}

names(listOfFiles2) <- files2

tf <- read.table('/dorothea/true_positives_ensembl.txt', header=F)

num <- vector()

for (i in 1:length(listOfFiles2)){
	temp <- listOfFiles2[[i]][((listOfFiles2[[i]]$V1 %in% tf$V1) | (listOfFiles2[[i]]$V2 %in% tf$V1)) & ((listOfFiles2[[i]]$V1 %in% tf$V2) | (listOfFiles2[[i]]$V2 %in% tf$V2)) ,]
	num[i] <- nrow(temp)/nrow(listOfFiles2[[i]])
}

names(num) <- names(listOfFiles2)

dat <- as.data.frame(num) %>% rownames_to_column("V2")

dat$Tissue[grepl("adipose",dat$V2)] <- "Adipose"
dat$Tissue[grepl("blood",dat$V2)] <- "Whole blood"
dat$Tissue[grepl("cmc",dat$V2)] <- "Prefrontal cortex"
dat$Tissue[grepl("common",dat$V2)] <- "Prefrontal cortex"
dat$Tissue[grepl("heart",dat$V2)] <- "Heart-left ventricle"
dat$Tissue[grepl("muscle",dat$V2)] <- "Skeletal muscle"
dat$Tissue[grepl("intestine",dat$V2)] <- "Small intestine"
dat$Tissue[grepl("spleen",dat$V2)] <- "Spleen"

dat$Adjustment[grepl("quant",dat$V2)] <- "None"
dat$Adjustment[grepl("ruv",dat$V2)] <- "RUVCorr"
dat$Adjustment[grepl("pc",dat$V2)] <- "PC"
dat$Adjustment[grepl("peer",dat$V2)] <- "PEER"
dat$Adjustment[grepl("known",dat$V2)] <- "Known\nCovariates"
dat$Adjustment[grepl("basic",dat$V2)] <- "None"
dat$Adjustment[grepl("confeti",dat$V2)] <- "CONFETI"


dat$color[dat$Tissue=="Prefrontal cortex"] <- "#FF6666"
dat$color[dat$Tissue=="Skeletal muscle"] <- "#003333"
dat$color[dat$Tissue=="Small intestine"] <- "#669900"
dat$color[dat$Tissue=="Adipose"] <- "#9999FF"
dat$color[dat$Tissue=="Spleen"] <- "#CCCC00"
dat$color[dat$Tissue=="Heart-left ventricle"] <- "#FF3399"
dat$color[dat$Tissue=="Whole blood"] <- "#3366CC"

dat3 <- dat
dat3$Cutoff <- "0.3"
############################################

setwd("/scripts/edge_list/cutoff_0.4")

files2 <- list.files(pattern = "*edge.list", full.names=FALSE)
#files2 <- files2[grepl("intestine",files2)]

listOfFiles2 <- lapply(files2, function(x) fread(x,drop=1)) 

for (i in 1:length(listOfFiles2)){
	listOfFiles2[[i]]$V1 <- gsub("\\..*","",listOfFiles2[[i]]$V1)
	listOfFiles2[[i]]$V2 <- gsub("\\..*","",listOfFiles2[[i]]$V2)
}

names(listOfFiles2) <- files2

tf <- read.table('/dorothea/true_positives_ensembl.txt', header=F)

num <- vector()

for (i in 1:length(listOfFiles2)){
	temp <- listOfFiles2[[i]][((listOfFiles2[[i]]$V1 %in% tf$V1) | (listOfFiles2[[i]]$V2 %in% tf$V1)) & ((listOfFiles2[[i]]$V1 %in% tf$V2) | (listOfFiles2[[i]]$V2 %in% tf$V2)) ,]
	num[i] <- nrow(temp)/nrow(listOfFiles2[[i]])
}

names(num) <- names(listOfFiles2)

dat <- as.data.frame(num) %>% rownames_to_column("V2")

dat$Tissue[grepl("adipose",dat$V2)] <- "Adipose"
dat$Tissue[grepl("blood",dat$V2)] <- "Whole blood"
dat$Tissue[grepl("cmc",dat$V2)] <- "Prefrontal cortex"
dat$Tissue[grepl("common",dat$V2)] <- "Prefrontal cortex"
dat$Tissue[grepl("heart",dat$V2)] <- "Heart-left ventricle"
dat$Tissue[grepl("muscle",dat$V2)] <- "Skeletal muscle"
dat$Tissue[grepl("intestine",dat$V2)] <- "Small intestine"
dat$Tissue[grepl("spleen",dat$V2)] <- "Spleen"

dat$Adjustment[grepl("quant",dat$V2)] <- "None"
dat$Adjustment[grepl("ruv",dat$V2)] <- "RUVCorr"
dat$Adjustment[grepl("pc",dat$V2)] <- "PC"
dat$Adjustment[grepl("peer",dat$V2)] <- "PEER"
dat$Adjustment[grepl("known",dat$V2)] <- "Known\nCovariates"
dat$Adjustment[grepl("basic",dat$V2)] <- "None"
dat$Adjustment[grepl("confeti",dat$V2)] <- "CONFETI"


dat$color[dat$Tissue=="Prefrontal cortex"] <- "#FF6666"
dat$color[dat$Tissue=="Skeletal muscle"] <- "#003333"
dat$color[dat$Tissue=="Small intestine"] <- "#669900"
dat$color[dat$Tissue=="Adipose"] <- "#9999FF"
dat$color[dat$Tissue=="Spleen"] <- "#CCCC00"
dat$color[dat$Tissue=="Heart-left ventricle"] <- "#FF3399"
dat$color[dat$Tissue=="Whole blood"] <- "#3366CC"

dat4 <- dat
dat4$Cutoff <- "0.4"
############################################


setwd("/scripts/edge_list/")

files2 <- list.files(pattern = "*edge.list", full.names=FALSE)
#files2 <- files2[grepl("intestine",files2)]

listOfFiles2 <- lapply(files2, function(x) fread(x,drop=1)) 

for (i in 1:length(listOfFiles2)){
	listOfFiles2[[i]]$V1 <- gsub("\\..*","",listOfFiles2[[i]]$V1)
	listOfFiles2[[i]]$V2 <- gsub("\\..*","",listOfFiles2[[i]]$V2)
}

names(listOfFiles2) <- files2

tf <- read.table('/dorothea/true_positives_ensembl.txt', header=F)

num <- vector()

for (i in 1:length(listOfFiles2)){
	temp <- listOfFiles2[[i]][((listOfFiles2[[i]]$V1 %in% tf$V1) | (listOfFiles2[[i]]$V2 %in% tf$V1)) & ((listOfFiles2[[i]]$V1 %in% tf$V2) | (listOfFiles2[[i]]$V2 %in% tf$V2)) ,]
	num[i] <- nrow(temp)/nrow(listOfFiles2[[i]])
}

names(num) <- names(listOfFiles2)

dat <- as.data.frame(num) %>% rownames_to_column("V2")

dat$Tissue[grepl("adipose",dat$V2)] <- "Adipose"
dat$Tissue[grepl("blood",dat$V2)] <- "Whole blood"
dat$Tissue[grepl("cmc",dat$V2)] <- "Prefrontal cortex"
dat$Tissue[grepl("common",dat$V2)] <- "Prefrontal cortex"
dat$Tissue[grepl("heart",dat$V2)] <- "Heart-left ventricle"
dat$Tissue[grepl("muscle",dat$V2)] <- "Skeletal muscle"
dat$Tissue[grepl("intestine",dat$V2)] <- "Small intestine"
dat$Tissue[grepl("spleen",dat$V2)] <- "Spleen"

dat$Adjustment[grepl("quant",dat$V2)] <- "None"
dat$Adjustment[grepl("ruv",dat$V2)] <- "RUVCorr"
dat$Adjustment[grepl("pc",dat$V2)] <- "PC"
dat$Adjustment[grepl("peer",dat$V2)] <- "PEER"
dat$Adjustment[grepl("known",dat$V2)] <- "Known\nCovariates"
dat$Adjustment[grepl("basic",dat$V2)] <- "None"
dat$Adjustment[grepl("confeti",dat$V2)] <- "CONFETI"


dat$color[dat$Tissue=="Prefrontal cortex"] <- "#FF6666"
dat$color[dat$Tissue=="Skeletal muscle"] <- "#003333"
dat$color[dat$Tissue=="Small intestine"] <- "#669900"
dat$color[dat$Tissue=="Adipose"] <- "#9999FF"
dat$color[dat$Tissue=="Spleen"] <- "#CCCC00"
dat$color[dat$Tissue=="Heart-left ventricle"] <- "#FF3399"
dat$color[dat$Tissue=="Whole blood"] <- "#3366CC"

dat5 <- dat
dat5$Cutoff <- "0.5"
############################################


setwd("/edge_list/cutoff_0.6")

files2 <- list.files(pattern = "*edge.list", full.names=FALSE)
#files2 <- files2[grepl("intestine",files2)]

listOfFiles2 <- lapply(files2, function(x) fread(x,drop=1)) 

for (i in 1:length(listOfFiles2)){
	listOfFiles2[[i]]$V1 <- gsub("\\..*","",listOfFiles2[[i]]$V1)
	listOfFiles2[[i]]$V2 <- gsub("\\..*","",listOfFiles2[[i]]$V2)
}

names(listOfFiles2) <- files2

tf <- read.table('/dorothea/true_positives_ensembl.txt', header=F)

num <- vector()

for (i in 1:length(listOfFiles2)){
	temp <- listOfFiles2[[i]][((listOfFiles2[[i]]$V1 %in% tf$V1) | (listOfFiles2[[i]]$V2 %in% tf$V1)) & ((listOfFiles2[[i]]$V1 %in% tf$V2) | (listOfFiles2[[i]]$V2 %in% tf$V2)) ,]
	num[i] <- nrow(temp)/nrow(listOfFiles2[[i]])
}

names(num) <- names(listOfFiles2)

dat <- as.data.frame(num) %>% rownames_to_column("V2")

dat$Tissue[grepl("adipose",dat$V2)] <- "Adipose"
dat$Tissue[grepl("blood",dat$V2)] <- "Whole blood"
dat$Tissue[grepl("cmc",dat$V2)] <- "Prefrontal cortex"
dat$Tissue[grepl("common",dat$V2)] <- "Prefrontal cortex"
dat$Tissue[grepl("heart",dat$V2)] <- "Heart-left ventricle"
dat$Tissue[grepl("muscle",dat$V2)] <- "Skeletal muscle"
dat$Tissue[grepl("intestine",dat$V2)] <- "Small intestine"
dat$Tissue[grepl("spleen",dat$V2)] <- "Spleen"

dat$Adjustment[grepl("quant",dat$V2)] <- "None"
dat$Adjustment[grepl("ruv",dat$V2)] <- "RUVCorr"
dat$Adjustment[grepl("pc",dat$V2)] <- "PC"
dat$Adjustment[grepl("peer",dat$V2)] <- "PEER"
dat$Adjustment[grepl("known",dat$V2)] <- "Known\nCovariates"
dat$Adjustment[grepl("basic",dat$V2)] <- "None"
dat$Adjustment[grepl("confeti",dat$V2)] <- "CONFETI"


dat$color[dat$Tissue=="Prefrontal cortex"] <- "#FF6666"
dat$color[dat$Tissue=="Skeletal muscle"] <- "#003333"
dat$color[dat$Tissue=="Small intestine"] <- "#669900"
dat$color[dat$Tissue=="Adipose"] <- "#9999FF"
dat$color[dat$Tissue=="Spleen"] <- "#CCCC00"
dat$color[dat$Tissue=="Heart-left ventricle"] <- "#FF3399"
dat$color[dat$Tissue=="Whole blood"] <- "#3366CC"

dat6 <- dat
dat6$Cutoff <- "0.6"
############################################



setwd("/scripts/edge_list/cutoff_0.7")

files2 <- list.files(pattern = "*edge.list", full.names=FALSE)
#files2 <- files2[grepl("intestine",files2)]

listOfFiles2 <- lapply(files2, function(x) fread(x,drop=1)) 

for (i in 1:length(listOfFiles2)){
	listOfFiles2[[i]]$V1 <- gsub("\\..*","",listOfFiles2[[i]]$V1)
	listOfFiles2[[i]]$V2 <- gsub("\\..*","",listOfFiles2[[i]]$V2)
}

names(listOfFiles2) <- files2

tf <- read.table('/dorothea/true_positives_ensembl.txt', header=F)

num <- vector()

for (i in 1:length(listOfFiles2)){
	temp <- listOfFiles2[[i]][((listOfFiles2[[i]]$V1 %in% tf$V1) | (listOfFiles2[[i]]$V2 %in% tf$V1)) & ((listOfFiles2[[i]]$V1 %in% tf$V2) | (listOfFiles2[[i]]$V2 %in% tf$V2)) ,]
	num[i] <- nrow(temp)/nrow(listOfFiles2[[i]])
}

names(num) <- names(listOfFiles2)

dat <- as.data.frame(num) %>% rownames_to_column("V2")

dat$Tissue[grepl("adipose",dat$V2)] <- "Adipose"
dat$Tissue[grepl("blood",dat$V2)] <- "Whole blood"
dat$Tissue[grepl("cmc",dat$V2)] <- "Prefrontal cortex"
dat$Tissue[grepl("common",dat$V2)] <- "Prefrontal cortex"
dat$Tissue[grepl("heart",dat$V2)] <- "Heart-left ventricle"
dat$Tissue[grepl("muscle",dat$V2)] <- "Skeletal muscle"
dat$Tissue[grepl("intestine",dat$V2)] <- "Small intestine"
dat$Tissue[grepl("spleen",dat$V2)] <- "Spleen"

dat$Adjustment[grepl("quant",dat$V2)] <- "None"
dat$Adjustment[grepl("ruv",dat$V2)] <- "RUVCorr"
dat$Adjustment[grepl("pc",dat$V2)] <- "PC"
dat$Adjustment[grepl("peer",dat$V2)] <- "PEER"
dat$Adjustment[grepl("known",dat$V2)] <- "Known\nCovariates"
dat$Adjustment[grepl("basic",dat$V2)] <- "None"
dat$Adjustment[grepl("confeti",dat$V2)] <- "CONFETI"


dat$color[dat$Tissue=="Prefrontal cortex"] <- "#FF6666"
dat$color[dat$Tissue=="Skeletal muscle"] <- "#003333"
dat$color[dat$Tissue=="Small intestine"] <- "#669900"
dat$color[dat$Tissue=="Adipose"] <- "#9999FF"
dat$color[dat$Tissue=="Spleen"] <- "#CCCC00"
dat$color[dat$Tissue=="Heart-left ventricle"] <- "#FF3399"
dat$color[dat$Tissue=="Whole blood"] <- "#3366CC"

dat7 <- dat
dat7$Cutoff <- "0.7"
############################################


master <- bind_rows(dat3,dat4,dat5,dat6,dat7)



pdf('~/www/percent_tf_edges_dorothea_all_cutoffs.pdf',width=10,height=11)
par(mfrow = c(1,1),mar=c(20,15,5,5))
boxplot(num ~ Adjustment + Cutoff, data=master, horizontal=FALSE,outline=FALSE,las=2,cex.lab=2, cex.axis=1.8, cex.main=2, cex.sub=2,boxcol="black", boxlty=1,whisklty=1,whisklwd=3,boxlwd=2,whiskcol="black",
	names=rep(c("CONFETI","Known\nCovariates","None","PC","PEER","RUVCorr"),5))
grid(col="grey",lty = "solid")
title(ylab="Proportion of edges",line=5,cex.lab=2)
title(xlab="Adjustment",line=10,cex.lab=2)
beeswarm(num ~ Adjustment + Cutoff, data=master, method="center",vertical=TRUE,pch = 16,las=2,pwcol=color,add=TRUE, cex=2)
legend("right",legend=unique(master$Tissue),col=unique(master$color),title="Tissue")
dev.off()


