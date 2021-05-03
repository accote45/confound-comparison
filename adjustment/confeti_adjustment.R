
library(confeti)
library(data.table)
library(tidyverse)
library(lrgpr)
library(Hmisc)
library(corrplot)
library(RColorBrewer)
library(variancePartition)
library(gplots)
library(biomaRt)
library(factoextra)
library(readxl)
library(gsubfn)

geno <- args[1]
expr.count <- args[2]
index <- rownames(expr.count)
geno2 <- geno[,colnames(geno) %in% index]
rm(geno)

expr2 <- expr.count[,order(colnames(expr.count))]
expr2[] <- lapply(expr2, function(x) as.numeric(as.character(x)))
geno3<- geno2[,order(colnames(geno2))]
rm(geno2)
geno3[] <- lapply(geno3, function(x) as.numeric(as.character(x)))

results = confeti(expr2,geno3)

# peaks are defined as gene contributions that are larger than 2 standard deviations
# save data table of peak genes per genetic IC factor
for (i in results$genetic_factors){
	dat <- as.data.frame(results$peaks[[i]])
	write.table(dat, paste0("gene_peaks_",i,".txt"), quote=FALSE)
}


ic_coefficients <- as.data.frame(results$A)
confounding <- ic_coefficients[c(results$confounding_factors),]
index <- names(confounding)

write.table(confounding,"ic_confounding_coefficients.txt", quote=FALSE)



results2=confeti(expr2,geno3,return_all=FALSE)
rm(geno3)

covs < args[3]
covs <- as.data.frame(t(covs)) %>% rownames_to_column("ID")
covs$SUBJID <- gsub("\\.","\\-",covs$ID)

index <- names(as.data.frame(results2$Kmx))
demos.fin <- subset(covs, covs$SUBJID %in% index)
demos.fin <- as.matrix(demos.fin[order(demos.fin$SUBJID),])

EV.1 = as.numeric(as.character(demos.fin[,2]))
EV.2 = as.numeric(as.character(demos.fin[,3]))
EV.3 = as.numeric(as.character(demos.fin[,4]))
EV.4 = as.numeric(as.character(demos.fin[,5]))
EV.5 = as.numeric(as.character(demos.fin[,6]))

expr.fin <- as.matrix(expr2)

residualsR<- list()
for (i in 1:nrow(expr.fin)){
	fit <- lrgpr( expr.fin[i,] ~ EV.1 + EV.2 + EV.3 + EV.4 + EV.5, decomp = svd(results2$Kmx))
	residualsR[[i]] <- residuals.lrgpr(fit)
	residuals.dat<- do.call(rbind, residualsR)
}

rownames(residuals.dat) <- rownames(expr.fin)
write.table(residuals.dat, "gtex_residual.gene_expression_confetiadj_wAnc.txt")
