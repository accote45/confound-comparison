



# gprofiler enrich, WGCNA

library(tidyverse)
library(plyr)
library(biomaRt)
library(purrr)
library(data.table)
library(biomaRt)
library(httr)
library(jsonlite)
library(gprofiler2)


args<-commandArgs(trailingOnly=TRUE)
mod=10
net="unsigned"

expr <- fread("expression.txt",stringsAsFactors=FALSE) %>% column_to_rownames(var = "SampleID")
expr <- t(expr)
rownames(expr) <- gsub("\\..*","",rownames(expr))

background_genes <- rownames(expr)

# read in modules for one method

dat <- list()
temp <- list()
listfiles <- list.files(pattern="module*")
listfiles <- listfiles[(!listfiles %in% paste0(net,"_pearson_min_mod_size",mod,"module_grey"))]
listfiles <- listfiles[-grep("unsigned_pearson_min_mod_size10module_eigengenes_10min.txt",listfiles)]
listfiles <- listfiles[-grep("unsigned_pearson_min_mod_size10module_assigned_master.txt",listfiles)]
listfiles <- listfiles[-grep("unsigned_pearson_min_mod_size10numbergenes_permodule.txt",listfiles)]
listfiles <- listfiles[-grep("unsigned_pearson_min_mod_size10wgcna_hclust_dendrogram_10min_mergedmodules.png",listfiles)]

#listfiles <- listfiles[-grep("unsigned_pearson_min_mod_size10wgcna_hclust_dendrogram_cmc1_10min_mergedmodules.png",listfiles)]


# replace ENS IDs with gene names
for (k in 1:length(listfiles)){
    temp[[k]]<- read.table(listfiles[k])
    dat[[k]] <- gsub("\\..*","",temp[[k]]$V1)
}

# query gprofiler -- 
enriched <- list()

for (k in 1:length(dat)){
temp <- gost(dat[[k]], organism = "hsapiens", ordered_query = FALSE, significant = FALSE, user_threshold = 1, correction_method = "fdr", domain_scope = "custom", custom_bg = background_genes, sources = c("GO", "KEGG", "REAC"))
enriched[[k]] <- temp$result
enriched[[k]]$term_name <- gsub(" ", "_", enriched[[k]]$term_name)
enriched[[k]]$parents <- gsub(" ", "_", enriched[[k]]$parents)
}

# save results files for each module
for (k in 1:length(enriched)){
    write.table(enriched[[k]], paste0("/GO_KEGG_REAC_enrich_results_", listfiles[[k]]), quote=FALSE, row.names=FALSE)
}


# count number of significant pathways for each module (per database)

num_enriched_mods <- list()

data <- c("GO", "KEGG", "REAC")
pvals <- as.numeric(c("22971", "337", "2250"))

for (k in 1:length(enriched)){
    num_enriched_mods[[k]] <- NaN*seq(3)
    for (i in 1:length(data))
    num_enriched_mods[[k]][[i]] <- sum(enriched[[k]][enriched[[k]]$source %like% data[[i]],]$p_value < (0.05/(length(listfiles)*pvals[[i]])), na.rm=TRUE)
}


dat <- as.data.frame(do.call(rbind, num_enriched_mods))
colnames(dat) <- data
rownames(dat) <- listfiles
write.table(dat, paste0("unsigned_pearson_min_mod_size10_sig_pathways_per_mod.txt"), quote=FALSE, row.names=TRUE)





