


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

coexpr="heart_known"

mods.raw <- read.table("FANTOM5_individual_networks/heart_adult.txt", stringsAsFactors=FALSE)

mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id_version", "hgnc_symbol", "gene_biotype"),
                       filters = "hgnc_symbol", values = unique(c(mods.raw$V1,mods.raw$V2)),
                       mart = mart, uniqueRows = TRUE)
new <- mods.raw 
new[] <- lapply(mods.raw[,1:2], function(x) Ensemble2HGNC$ensembl_gene_id_version[match(x, Ensemble2HGNC$hgnc_symbol)])
new$V3 <- mods.raw$V3
new$V1 <- gsub("\\..*","",new$V1)
new$V2 <- gsub("\\..*","",new$V2)

setwd(paste0("/modules/",coexpr))

listfiles <- list.files(pattern=glob2rx("unsigned_pearson_min_mod_size10module*"))
listfiles <- listfiles[!(listfiles %in% c("unsigned_pearson_min_mod_size10module_assigned_master.txt","unsigned_pearson_min_mod_size10module_eigengenes_10min.txt"))]

# read in coexpression modules
coexpr_mods <- list()
for (k in 1:length(listfiles)){
	coexpr_mods[[k]] <- read.table(text=readLines(listfiles[k]), header=FALSE, quote="")
	coexpr_mods[[k]] <- gsub("\\..*","",coexpr_mods[[k]]$V1) 	
}
names(coexpr_mods) <- listfiles

sizes <- read.table("gene_count_per_file.txt")
domain_size <- sizes[sizes$V1==coexpr,]$V2

# get intersection of genes for each module and each regulator group across cutoffs

cutoffs <- c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)

for (cutoff in cutoffs) {
	ref <- new[new$V3>cutoff,]
	assign(paste0("mods.",cutoff),split(ref$V2,ref$V1))
	gom.obj <- newGOM(eval(parse(text=paste0("mods.",cutoff))), coexpr_mods, genome.size=domain_size)
	pvals <- getMatrix(gom.obj, name="pval") %>% as.data.frame %>% rownames_to_column("TF")
	OR <- getMatrix(gom.obj, "odds.ratio") %>% as.data.frame %>% rownames_to_column("TF")
	pvalsfin <- pvals %>% gather(module,p_values,colnames(pvals)[-1])
	ORfin <- OR %>% gather(module,OR,colnames(OR)[-1])
	assign(paste0("master.",cutoff),merge(ORfin,pvalsfin,by=c("TF","module")))
}


d1.aucodds = function(i){
	files[[i]]$q_value<- p.adjust(files[[i]]$p_values,method="holm")
	master.sig <- files[[i]][files[[i]]$q_value<0.1,]
	master.sig <- as.data.table(master.sig)

	master.dat <- master.sig[master.sig[, .I[OR==max(OR)], by=TF]$V1]
	# calculate percentage of regulator sets with OR > cut off and AUC

	x <- seq(0, 3, length.out=1000)
		d1 <- data.frame(proportion=rep(0,length(x)), log10_OR=rep(0,length(x)))
		for (z in 1:length(x)){
			temp <- log10(master.dat$OR)>=x[z]
			d1[z,] <- c(length(temp[temp=="TRUE"])/length(get(paste0("mods.",names(files)[[i]]))), x[z])
		}
	d1.aucodds <- trapz(d1$log10_OR, d1$proportion)
	return(d1.aucodds)
}

d1.aucodds.dat <- lapply(seq_along(files), d1.aucodds)

