

library(tidyverse)
library(RUVcorr)
library(data.table)


ge.mat <- args[1]


library(tidyverse)
library(RUVcorr)
library(data.table)



nc_index <- empNegativeControls(as.matrix(ge.mat), exclude=1, nc=2000)

pdf(IQR_vs_mean_nc_genes.pdf)
genePlot(ge.mat, index=nc_index, legend="Negative Control Genes", title="IQR-Mean Plot")
dev.off()


library(snowfall)

k <- c(1,2,3,4,5)
nu <- c(0,500,1000,5000)
k.nu.matrix <- cbind(rep(k, each=4), rep(nu, 5))
k.nu.matrix <- as.list(as.data.frame(t(k.nu.matrix)))
sfInit(parallel=TRUE, cpus=4)

sfLibrary(RUVcorr)

sfExport("ge.mat", "k.nu.matrix", "nc_index")
expr_AllRUV <- sfLapply(k.nu.matrix, function(x)
RUVNaiveRidge(ge.mat, center=TRUE, nc_index, x[2], x[1]))
sfStop()

# sodium channel genes
na_genes <- c("ENSG00000144285", "ENSG00000153253", "ENSG00000007314", "ENSG00000183873", "ENSG00000136546", "ENSG00000196876", "ENSG00000168356",
"ENSG00000105711", "ENSG00000149575", "ENSG00000166257", "ENSG00000177098")
na_index <- which(is.element(colnames(ge.mat), na_genes))

cor_AllRUV_na <- lapply(expr_AllRUV, function(x) cor(x[,na_index]))
cor_Raw_na <- cor(ge.mat[,na_index])

pdf("sodium_channel_genes_param_test.pdf")
par(mfrow=c(2,2))
lapply(1:4, function(i) histogramPlot(cor_AllRUV_na[seq(0,18,4)+i], cor_Raw_na,
title=paste("nu=", nu[i]),
legend=c(paste("k=", k), "Raw")))
dev.off()


# MHC genes
na_genes <- c("ENSG00000206503", "ENSG00000234745", "ENSG00000204525", "ENSG00000204257", "ENSG00000242574", "ENSG00000204252", "ENSG00000241106","ENSG00000231389","ENSG00000231461",
	"ENSG00000237398","ENSG00000223865","ENSG00000224557","ENSG00000196735","ENSG00000237541","ENSG00000179344","ENSG00000232629", "ENSG00000226030", "ENSG00000204287", "ENSG00000196126",
	"ENSG00000227442","ENSG00000196101","ENSG00000227357","ENSG00000198502","ENSG00000229391","ENSG00000227099","ENSG00000233697","ENSG00000196301","ENSG00000204592","ENSG00000204642",
	"ENSG00000204632","ENSG00000206341","ENSG00000204622","ENSG00000230795","ENSG00000243753","ENSG00000224372","ENSG00000261548","ENSG00000225851","ENSG00000231130","ENSG00000228078",
	"ENSG00000181126","ENSG00000235290","ENSG00000235301")
na_index <- which(is.element(colnames(ge.mat), na_genes))

cor_AllRUV_na <- lapply(expr_AllRUV, function(x) cor(x[,na_index]))
cor_Raw_na <- cor(ge.mat[,na_index])

pdf("mhc_genes_param_test.pdf")
par(mfrow=c(2,2))
lapply(1:4, function(i) histogramPlot(cor_AllRUV_na[seq(0,18,4)+i], cor_Raw_na,
title=paste("nu=", nu[i]),
legend=c(paste("k=", k), "Raw")))
dev.off()

library(biomaRt)
rib_genes <- read.table("rib_genes.csv", header=TRUE)
# Define biomart object
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Query biomart
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                       filters = "hgnc_symbol", values = rib_genes$gene,
                       mart = mart, uniqueRows = TRUE)

na_genes<-Ensemble2HGNC$ensembl_gene_id
na_index <- which(is.element(colnames(ge.mat), na_genes))

cor_AllRUV_na <- lapply(expr_AllRUV, function(x) cor(x[,na_index]))
cor_Raw_na <- cor(ge.mat[,na_index])

pdf("rib_genes_param_test.pdf")
par(mfrow=c(2,2))
lapply(1:4, function(i) histogramPlot(cor_AllRUV_na[seq(0,18,4)+i], cor_Raw_na,
title=paste("nu=", nu[i]),
legend=c(paste("k=", k), "Raw")))
dev.off()


bg_index <- background(ge.mat, nBG=1000, exclude=na_index, nc_index=nc_index)
cor_AllRUV_bg <- lapply(expr_AllRUV, function(x) cor(x[,bg_index]))
cor_Raw_bg <- cor(ge.mat[,bg_index])
pdf("random_genes_param_test.pdf")
par(mfrow=c(2,2))
lapply(1:4, function(i) histogramPlot(cor_AllRUV_bg[seq(0,18,4)+i], cor_Raw_bg,
title=paste("nu=", nu[i]),
legend=c(paste("k=", k), "Raw")))
dev.off()

pdf("k4_boxplot_genes_param_test.pdf"))
par(mfrow=c(2,2))
lapply(1:4, function(i) RLEPlot(ge.mat, expr_AllRUV[[12+i]],
name=c("Raw", "RUV"), title=paste("nu=", nu[i]),
method="IQR.boxplots"))
dev.off()
