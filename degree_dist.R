



library(tidyverse)
library(plyr)
library(data.table)
library(purrr)
library(igraph)
library(reshape)

# read in module files
args<-commandArgs(trailingOnly=TRUE)
folder=args[1]
file=args[2]


dat <- fread(paste0("expr/",folder,"/",file)) %>% column_to_rownames("V1")

graph.obj<- graph_from_adjacency_matrix(as.matrix(cor(t(dat), method="pearson")), mode="undirected", weighted=TRUE, diag=FALSE)
graph.obj <- simplify(graph.obj, remove.multiple=TRUE, remove.loops=TRUE)

# Remove edges below absolute Pearson correlation 0.5
graph.obj <- delete_edges(graph.obj, E(graph.obj)[which(abs(E(graph.obj)$weight)<0.5)])
# Remove any vertices remaining that have no edges
graph.obj <- delete.vertices(graph.obj, degree(graph.obj)==0)

# calculate number of edges per gene
degree.cum.basic <- degree_distribution(graph.obj,cumulative=TRUE)

write.table(degree.cum.basic, paste0(file,"_degree.cum.basic"), quote=FALSE)
