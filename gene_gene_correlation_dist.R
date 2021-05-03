

library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)

basic <- args[1]
known <- args[2]
pc <- args[3]
peer <- args[4]
ruvcorr <- args[5]
confeti <- args[6]

basic.fin <- basic[sample(nrow(basic),5000),]
basic2 <- cor(t(as.matrix(basic.fin)))
basic3 <- melt(basic2)
#basic4 <- basic3[!duplicated(basic3)]
rm(basic2)
rm(basic.fin)
rm(basic)


known.fin <- known[sample(nrow(known),5000),]
known2 <- cor(t(as.matrix(known.fin)))
known3 <- melt(known2)
rm(known2)
rm(known.fin)
rm(known)


pc.fin <- t(pc)[sample(nrow(t(pc)),5000),]
pc2 <- cor(t(as.matrix(pc.fin)))
pc3 <- melt(pc2)
rm(pc2)
rm(pc.fin)
rm(pc)


peer.fin <- t(peer)[sample(nrow(t(peer)),5000),]
peer2 <- cor(t(as.matrix(peer.fin)))
peer3 <- melt(peer2)
rm(peer2)
rm(peer.fin)
rm(peer)


ruvcorr.fin <- ruvcorr[sample(nrow(ruvcorr),5000),]
ruvcorr2 <- cor(t(as.matrix(ruvcorr.fin)))
ruvcorr3 <- melt(ruvcorr2)
rm(ruvcorr2)
rm(ruvcorr.fin)
rm(ruvcorr)


confeti.fin <- confeti[sample(nrow(confeti),5000),]
confeti2 <- cor(t(as.matrix(confeti.fin)))
confeti3 <- melt(confeti2)
rm(confeti2)
rm(confeti.fin)
rm(confeti)


# plot distribution of gene-gene correlations across whole genome

master <- data.frame(basic3$value,known3$value,pc3$value,peer3$value,ruvcorr3$value,confeti3$value)
colnames(master) <- c("None","Known Covariates","PC","PEER","RUVCorr","CONFETI")

master.fin <- melt(master)

pdf("random_gene_corr.pdf")
p <- master.fin %>%
  ggplot( aes(x=value, color=variable)) +
    geom_density() + theme_bw()+ labs(color="Adjustment Method") + xlab("Pearson correlation coefficient") + ylab("Density") + scale_color_manual(values=c("Known Covariates"="#6699CC","None"="#FF99CC","RUVCorr"="#CC3300","CONFETI"="#666600","PC"="#9900FF","PEER"="#00CCCC")) 
p
dev.off()
