


# overlap between modules after each adjustment method 


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
library(gplots)
library(RColorBrewer)
library(patchwork)

tissue="small_intestine"

folders <- c(paste0(tissue,"_ruvcorr"),paste0(tissue,"_pc"),paste0(tissue,"_known"),paste0(tissue,"_basic"),paste0(tissue,"_confeti"),paste0(tissue,"_peer"))


for (folder in folders){
	setwd(paste0("/sc/arion/projects/psychgen2/cotea02/confound_comparisons/results/modules/",folder))
	listfiles <- list.files(pattern=glob2rx("unsigned_pearson_min_mod_size10module*"))
	listfiles <- listfiles[!(listfiles %in% c("unsigned_pearson_min_mod_size10module_assigned_master.txt","unsigned_pearson_min_mod_size10module_eigengenes_10min.txt"))]
	coexpr_mods <- list()
	for (k in 1:length(listfiles)){
		coexpr_mods[[k]] <- read.table(text=readLines(listfiles[k]), header=FALSE, quote="")
		coexpr_mods[[k]] <- gsub("\\..*","",coexpr_mods[[k]]$V1) 	
	}
	names(coexpr_mods) <- paste0(listfiles,"_",folder)
	assign(paste0("mods.",folder),coexpr_mods)
}


# get intersection of genes for every module pair across adjustment methods
master <- c(eval(parse(text=paste0("mods.",tissue,"_ruvcorr"))),eval(parse(text=paste0("mods.",tissue,"_pc"))),eval(parse(text=paste0("mods.",tissue,"_known"))),eval(parse(text=paste0("mods.",tissue,"_basic"))),eval(parse(text=paste0("mods.",tissue,"_confeti"))),eval(parse(text=paste0("mods.",tissue,"_peer"))))


gom.obj <- newGOM(master, master)
jaccard <- getMatrix(gom.obj, name="Jaccard") %>% as.data.frame %>% rownames_to_column("module")

dat <- melt(jaccard) 
dat$ID1 <- gsub(".*\\_","",dat$module)
dat$ID2 <- gsub(".*\\_","",dat$variable)

master <- dat %>% group_by(ID1,ID2) %>% dplyr::summarise(n=sum(value>0.5))

# make within adj comparisons set to null
master[master$ID1==master$ID2,]$n <- NA

master$ID1[master$ID1=="basic"] <- "None"
master$ID1[master$ID1=="confeti"] <- "CONFETI"
master$ID1[master$ID1=="known"] <- "Known\nCovariates"
master$ID1[master$ID1=="pc"] <- "PC"
master$ID1[master$ID1=="peer"] <- "PEER"
master$ID1[master$ID1=="ruvcorr"] <- "RUVCorr"

master$ID2[master$ID2=="basic"] <- "None"
master$ID2[master$ID2=="confeti"] <- "CONFETI"
master$ID2[master$ID2=="known"] <- "Known\nCovariates"
master$ID2[master$ID2=="pc"] <- "PC"
master$ID2[master$ID2=="peer"] <- "PEER"
master$ID2[master$ID2=="ruvcorr"] <- "RUVCorr"

master.intestine <- master



####################################
tissue="spleen"

folders <- c(paste0(tissue,"_ruvcorr"),paste0(tissue,"_pc"),paste0(tissue,"_known"),paste0(tissue,"_basic"),paste0(tissue,"_confeti"),paste0(tissue,"_peer"))


for (folder in folders){
	setwd(paste0("/sc/arion/projects/psychgen2/cotea02/confound_comparisons/results/modules/",folder))
	listfiles <- list.files(pattern=glob2rx("unsigned_pearson_min_mod_size10module*"))
	listfiles <- listfiles[!(listfiles %in% c("unsigned_pearson_min_mod_size10module_assigned_master.txt","unsigned_pearson_min_mod_size10module_eigengenes_10min.txt"))]
	coexpr_mods <- list()
	for (k in 1:length(listfiles)){
		coexpr_mods[[k]] <- read.table(text=readLines(listfiles[k]), header=FALSE, quote="")
		coexpr_mods[[k]] <- gsub("\\..*","",coexpr_mods[[k]]$V1) 	
	}
	names(coexpr_mods) <- paste0(listfiles,"_",folder)
	assign(paste0("mods.",folder),coexpr_mods)
}


# get intersection of genes for every module pair across adjustment methods
master <- c(eval(parse(text=paste0("mods.",tissue,"_ruvcorr"))),eval(parse(text=paste0("mods.",tissue,"_pc"))),eval(parse(text=paste0("mods.",tissue,"_known"))),eval(parse(text=paste0("mods.",tissue,"_basic"))),eval(parse(text=paste0("mods.",tissue,"_confeti"))),eval(parse(text=paste0("mods.",tissue,"_peer"))))


gom.obj <- newGOM(master, master)
jaccard <- getMatrix(gom.obj, name="Jaccard") %>% as.data.frame %>% rownames_to_column("module")

dat <- melt(jaccard) 
dat$ID1 <- gsub(".*\\_","",dat$module)
dat$ID2 <- gsub(".*\\_","",dat$variable)

master <- dat %>% group_by(ID1,ID2) %>% dplyr::summarise(n=sum(value>0.5))


# make within adj comparisons set to null
master[master$ID1==master$ID2,]$n <- NA

master$ID1[master$ID1=="basic"] <- "None"
master$ID1[master$ID1=="confeti"] <- "CONFETI"
master$ID1[master$ID1=="known"] <- "Known\nCovariates"
master$ID1[master$ID1=="pc"] <- "PC"
master$ID1[master$ID1=="peer"] <- "PEER"
master$ID1[master$ID1=="ruvcorr"] <- "RUVCorr"

master$ID2[master$ID2=="basic"] <- "None"
master$ID2[master$ID2=="confeti"] <- "CONFETI"
master$ID2[master$ID2=="known"] <- "Known\nCovariates"
master$ID2[master$ID2=="pc"] <- "PC"
master$ID2[master$ID2=="peer"] <- "PEER"
master$ID2[master$ID2=="ruvcorr"] <- "RUVCorr"

master.spleen <- master



####################################
tissue="heart"

folders <- c(paste0(tissue,"_ruvcorr"),paste0(tissue,"_pc"),paste0(tissue,"_known"),paste0(tissue,"_basic"),paste0(tissue,"_confeti"),paste0(tissue,"_peer"))


for (folder in folders){
	setwd(paste0("/sc/arion/projects/psychgen2/cotea02/confound_comparisons/results/modules/",folder))
	listfiles <- list.files(pattern=glob2rx("unsigned_pearson_min_mod_size10module*"))
	listfiles <- listfiles[!(listfiles %in% c("unsigned_pearson_min_mod_size10module_assigned_master.txt","unsigned_pearson_min_mod_size10module_eigengenes_10min.txt"))]
	coexpr_mods <- list()
	for (k in 1:length(listfiles)){
		coexpr_mods[[k]] <- read.table(text=readLines(listfiles[k]), header=FALSE, quote="")
		coexpr_mods[[k]] <- gsub("\\..*","",coexpr_mods[[k]]$V1) 	
	}
	names(coexpr_mods) <- paste0(listfiles,"_",folder)
	assign(paste0("mods.",folder),coexpr_mods)
}


# get intersection of genes for every module pair across adjustment methods
master <- c(eval(parse(text=paste0("mods.",tissue,"_ruvcorr"))),eval(parse(text=paste0("mods.",tissue,"_pc"))),eval(parse(text=paste0("mods.",tissue,"_known"))),eval(parse(text=paste0("mods.",tissue,"_basic"))),eval(parse(text=paste0("mods.",tissue,"_confeti"))),eval(parse(text=paste0("mods.",tissue,"_peer"))))


gom.obj <- newGOM(master, master)
jaccard <- getMatrix(gom.obj, name="Jaccard") %>% as.data.frame %>% rownames_to_column("module")

dat <- melt(jaccard) 
dat$ID1 <- gsub(".*\\_","",dat$module)
dat$ID2 <- gsub(".*\\_","",dat$variable)

master <- dat %>% group_by(ID1,ID2) %>% dplyr::summarise(n=sum(value>0.5))


# make within adj comparisons set to null
master[master$ID1==master$ID2,]$n <- NA

master$ID1[master$ID1=="basic"] <- "None"
master$ID1[master$ID1=="confeti"] <- "CONFETI"
master$ID1[master$ID1=="known"] <- "Known\nCovariates"
master$ID1[master$ID1=="pc"] <- "PC"
master$ID1[master$ID1=="peer"] <- "PEER"
master$ID1[master$ID1=="ruvcorr"] <- "RUVCorr"

master$ID2[master$ID2=="basic"] <- "None"
master$ID2[master$ID2=="confeti"] <- "CONFETI"
master$ID2[master$ID2=="known"] <- "Known\nCovariates"
master$ID2[master$ID2=="pc"] <- "PC"
master$ID2[master$ID2=="peer"] <- "PEER"
master$ID2[master$ID2=="ruvcorr"] <- "RUVCorr"

master.heart <- master





####################################
tissue="adipose"

folders <- c(paste0(tissue,"_ruvcorr"),paste0(tissue,"_pc"),paste0(tissue,"_known"),paste0(tissue,"_basic"),paste0(tissue,"_confeti"),paste0(tissue,"_peer"))


for (folder in folders){
	setwd(paste0("/sc/arion/projects/psychgen2/cotea02/confound_comparisons/results/modules/",folder))
	listfiles <- list.files(pattern=glob2rx("unsigned_pearson_min_mod_size10module*"))
	listfiles <- listfiles[!(listfiles %in% c("unsigned_pearson_min_mod_size10module_assigned_master.txt","unsigned_pearson_min_mod_size10module_eigengenes_10min.txt"))]
	coexpr_mods <- list()
	for (k in 1:length(listfiles)){
		coexpr_mods[[k]] <- read.table(text=readLines(listfiles[k]), header=FALSE, quote="")
		coexpr_mods[[k]] <- gsub("\\..*","",coexpr_mods[[k]]$V1) 	
	}
	names(coexpr_mods) <- paste0(listfiles,"_",folder)
	assign(paste0("mods.",folder),coexpr_mods)
}


# get intersection of genes for every module pair across adjustment methods
master <- c(eval(parse(text=paste0("mods.",tissue,"_ruvcorr"))),eval(parse(text=paste0("mods.",tissue,"_pc"))),eval(parse(text=paste0("mods.",tissue,"_known"))),eval(parse(text=paste0("mods.",tissue,"_basic"))),eval(parse(text=paste0("mods.",tissue,"_confeti"))),eval(parse(text=paste0("mods.",tissue,"_peer"))))


gom.obj <- newGOM(master, master)
jaccard <- getMatrix(gom.obj, name="Jaccard") %>% as.data.frame %>% rownames_to_column("module")

dat <- melt(jaccard) 
dat$ID1 <- gsub(".*\\_","",dat$module)
dat$ID2 <- gsub(".*\\_","",dat$variable)

master <- dat %>% group_by(ID1,ID2) %>% dplyr::summarise(n=sum(value>0.5))

# make within adj comparisons set to null
master[master$ID1==master$ID2,]$n <- NA

master$ID1[master$ID1=="basic"] <- "None"
master$ID1[master$ID1=="confeti"] <- "CONFETI"
master$ID1[master$ID1=="known"] <- "Known\nCovariates"
master$ID1[master$ID1=="pc"] <- "PC"
master$ID1[master$ID1=="peer"] <- "PEER"
master$ID1[master$ID1=="ruvcorr"] <- "RUVCorr"

master$ID2[master$ID2=="basic"] <- "None"
master$ID2[master$ID2=="confeti"] <- "CONFETI"
master$ID2[master$ID2=="known"] <- "Known\nCovariates"
master$ID2[master$ID2=="pc"] <- "PC"
master$ID2[master$ID2=="peer"] <- "PEER"
master$ID2[master$ID2=="ruvcorr"] <- "RUVCorr"

master.adipose <- master




####################################
tissue="muscle"

folders <- c(paste0(tissue,"_ruvcorr"),paste0(tissue,"_pc"),paste0(tissue,"_known"),paste0(tissue,"_basic"),paste0(tissue,"_confeti"),paste0(tissue,"_peer"))


for (folder in folders){
	setwd(paste0("/sc/arion/projects/psychgen2/cotea02/confound_comparisons/results/modules/",folder))
	listfiles <- list.files(pattern=glob2rx("unsigned_pearson_min_mod_size10module*"))
	listfiles <- listfiles[!(listfiles %in% c("unsigned_pearson_min_mod_size10module_assigned_master.txt","unsigned_pearson_min_mod_size10module_eigengenes_10min.txt"))]
	coexpr_mods <- list()
	for (k in 1:length(listfiles)){
		coexpr_mods[[k]] <- read.table(text=readLines(listfiles[k]), header=FALSE, quote="")
		coexpr_mods[[k]] <- gsub("\\..*","",coexpr_mods[[k]]$V1) 	
	}
	names(coexpr_mods) <- paste0(listfiles,"_",folder)
	assign(paste0("mods.",folder),coexpr_mods)
}


# get intersection of genes for every module pair across adjustment methods
master <- c(eval(parse(text=paste0("mods.",tissue,"_ruvcorr"))),eval(parse(text=paste0("mods.",tissue,"_pc"))),eval(parse(text=paste0("mods.",tissue,"_known"))),eval(parse(text=paste0("mods.",tissue,"_basic"))),eval(parse(text=paste0("mods.",tissue,"_confeti"))),eval(parse(text=paste0("mods.",tissue,"_peer"))))


gom.obj <- newGOM(master, master)
jaccard <- getMatrix(gom.obj, name="Jaccard") %>% as.data.frame %>% rownames_to_column("module")

dat <- melt(jaccard) 
dat$ID1 <- gsub(".*\\_","",dat$module)
dat$ID2 <- gsub(".*\\_","",dat$variable)

master <- dat %>% group_by(ID1,ID2) %>% dplyr::summarise(n=sum(value>0.5))


# make within adj comparisons set to null
master[master$ID1==master$ID2,]$n <- NA

master$ID1[master$ID1=="basic"] <- "None"
master$ID1[master$ID1=="confeti"] <- "CONFETI"
master$ID1[master$ID1=="known"] <- "Known\nCovariates"
master$ID1[master$ID1=="pc"] <- "PC"
master$ID1[master$ID1=="peer"] <- "PEER"
master$ID1[master$ID1=="ruvcorr"] <- "RUVCorr"

master$ID2[master$ID2=="basic"] <- "None"
master$ID2[master$ID2=="confeti"] <- "CONFETI"
master$ID2[master$ID2=="known"] <- "Known\nCovariates"
master$ID2[master$ID2=="pc"] <- "PC"
master$ID2[master$ID2=="peer"] <- "PEER"
master$ID2[master$ID2=="ruvcorr"] <- "RUVCorr"

master.muscle <- master






####################################
tissue="blood"

folders <- c(paste0(tissue,"_ruvcorr"),paste0(tissue,"_pc"),paste0(tissue,"_known"),paste0(tissue,"_basic"),paste0(tissue,"_confeti"),paste0(tissue,"_peer"))


for (folder in folders){
	setwd(paste0("/sc/arion/projects/psychgen2/cotea02/confound_comparisons/results/modules/",folder))
	listfiles <- list.files(pattern=glob2rx("unsigned_pearson_min_mod_size10module*"))
	listfiles <- listfiles[!(listfiles %in% c("unsigned_pearson_min_mod_size10module_assigned_master.txt","unsigned_pearson_min_mod_size10module_eigengenes_10min.txt"))]
	coexpr_mods <- list()
	for (k in 1:length(listfiles)){
		coexpr_mods[[k]] <- read.table(text=readLines(listfiles[k]), header=FALSE, quote="")
		coexpr_mods[[k]] <- gsub("\\..*","",coexpr_mods[[k]]$V1) 	
	}
	names(coexpr_mods) <- paste0(listfiles,"_",folder)
	assign(paste0("mods.",folder),coexpr_mods)
}


# get intersection of genes for every module pair across adjustment methods
master <- c(eval(parse(text=paste0("mods.",tissue,"_ruvcorr"))),eval(parse(text=paste0("mods.",tissue,"_pc"))),eval(parse(text=paste0("mods.",tissue,"_known"))),eval(parse(text=paste0("mods.",tissue,"_basic"))),eval(parse(text=paste0("mods.",tissue,"_confeti"))),eval(parse(text=paste0("mods.",tissue,"_peer"))))


gom.obj <- newGOM(master, master)
jaccard <- getMatrix(gom.obj, name="Jaccard") %>% as.data.frame %>% rownames_to_column("module")

dat <- melt(jaccard) 
dat$ID1 <- gsub(".*\\_","",dat$module)
dat$ID2 <- gsub(".*\\_","",dat$variable)

master <- dat %>% group_by(ID1,ID2) %>% dplyr::summarise(n=sum(value>0.5))


# make within adj comparisons set to null
master[master$ID1==master$ID2,]$n <- NA

master$ID1[master$ID1=="basic"] <- "None"
master$ID1[master$ID1=="confeti"] <- "CONFETI"
master$ID1[master$ID1=="known"] <- "Known\nCovariates"
master$ID1[master$ID1=="pc"] <- "PC"
master$ID1[master$ID1=="peer"] <- "PEER"
master$ID1[master$ID1=="ruvcorr"] <- "RUVCorr"

master$ID2[master$ID2=="basic"] <- "None"
master$ID2[master$ID2=="confeti"] <- "CONFETI"
master$ID2[master$ID2=="known"] <- "Known\nCovariates"
master$ID2[master$ID2=="pc"] <- "PC"
master$ID2[master$ID2=="peer"] <- "PEER"
master$ID2[master$ID2=="ruvcorr"] <- "RUVCorr"

master.blood<- master















master.intestine$ID1[master.intestine$ID1=="None"] <- "None (25)"
master.intestine$ID1[master.intestine$ID1=="CONFETI"] <- "CONFETI (24)"
master.intestine$ID1[master.intestine$ID1=="Known\nCovariates"] <- "Known\nCovariates (35)"
master.intestine$ID1[master.intestine$ID1=="PC"] <- "PC (139)"
master.intestine$ID1[master.intestine$ID1=="PEER"] <- "PEER (117)"
master.intestine$ID1[master.intestine$ID1=="RUVCorr"] <- "RUVCorr (55)"

master.intestine$ID2[master.intestine$ID2=="None"] <- "None (25)"
master.intestine$ID2[master.intestine$ID2=="CONFETI"] <- "CONFETI (24)"
master.intestine$ID2[master.intestine$ID2=="Known\nCovariates"] <- "Known\nCovariates (35)"
master.intestine$ID2[master.intestine$ID2=="PC"] <- "PC (139)"
master.intestine$ID2[master.intestine$ID2=="PEER"] <- "PEER (117)"
master.intestine$ID2[master.intestine$ID2=="RUVCorr"] <- "RUVCorr (55)"





master.spleen$ID1[master.spleen$ID1=="None"] <- "None (56)"
master.spleen$ID1[master.spleen$ID1=="CONFETI"] <- "CONFETI (187)"
master.spleen$ID1[master.spleen$ID1=="Known\nCovariates"] <- "Known\nCovariates (90)"
master.spleen$ID1[master.spleen$ID1=="PC"] <- "PC (112)"
master.spleen$ID1[master.spleen$ID1=="PEER"] <- "PEER (133)"
master.spleen$ID1[master.spleen$ID1=="RUVCorr"] <- "RUVCorr (82)"

master.spleen$ID2[master.spleen$ID2=="None"] <- "None (56)"
master.spleen$ID2[master.spleen$ID2=="CONFETI"] <- "CONFETI (187)"
master.spleen$ID2[master.spleen$ID2=="Known\nCovariates"] <- "Known\nCovariates (90)"
master.spleen$ID2[master.spleen$ID2=="PC"] <- "PC (112)"
master.spleen$ID2[master.spleen$ID2=="PEER"] <- "PEER (133)"
master.spleen$ID2[master.spleen$ID2=="RUVCorr"] <- "RUVCorr (82)"




master.heart$ID1[master.heart$ID1=="None"] <- "None (33)"
master.heart$ID1[master.heart$ID1=="CONFETI"] <- "CONFETI (69)"
master.heart$ID1[master.heart$ID1=="Known\nCovariates"] <- "Known\nCovariates (57)"
master.heart$ID1[master.heart$ID1=="PC"] <- "PC (140)"
master.heart$ID1[master.heart$ID1=="PEER"] <- "PEER (89)"
master.heart$ID1[master.heart$ID1=="RUVCorr"] <- "RUVCorr (114)"

master.heart$ID2[master.heart$ID2=="None"] <- "None (33)"
master.heart$ID2[master.heart$ID2=="CONFETI"] <- "CONFETI (69)"
master.heart$ID2[master.heart$ID2=="Known\nCovariates"] <- "Known\nCovariates (57)"
master.heart$ID2[master.heart$ID2=="PC"] <- "PC (140)"
master.heart$ID2[master.heart$ID2=="PEER"] <- "PEER (89)"
master.heart$ID2[master.heart$ID2=="RUVCorr"] <- "RUVCorr (114)"





master.muscle$ID1[master.muscle$ID1=="None"] <- "None (74)"
master.muscle$ID1[master.muscle$ID1=="CONFETI"] <- "CONFETI (117)"
master.muscle$ID1[master.muscle$ID1=="Known\nCovariates"] <- "Known\nCovariates (120)"
master.muscle$ID1[master.muscle$ID1=="PC"] <- "PC (122)"
master.muscle$ID1[master.muscle$ID1=="PEER"] <- "PEER (113)"
master.muscle$ID1[master.muscle$ID1=="RUVCorr"] <- "RUVCorr (113)"

master.muscle$ID2[master.muscle$ID2=="None"] <- "None (74)"
master.muscle$ID2[master.muscle$ID2=="CONFETI"] <- "CONFETI (117)"
master.muscle$ID2[master.muscle$ID2=="Known\nCovariates"] <- "Known\nCovariates (120)"
master.muscle$ID2[master.muscle$ID2=="PC"] <- "PC (122)"
master.muscle$ID2[master.muscle$ID2=="PEER"] <- "PEER (113)"
master.muscle$ID2[master.muscle$ID2=="RUVCorr"] <- "RUVCorr (113)"




master.adipose$ID1[master.adipose$ID1=="None"] <- "None (79)"
master.adipose$ID1[master.adipose$ID1=="CONFETI"] <- "CONFETI (93)"
master.adipose$ID1[master.adipose$ID1=="Known\nCovariates"] <- "Known\nCovariates (102)"
master.adipose$ID1[master.adipose$ID1=="PC"] <- "PC (136)"
master.adipose$ID1[master.adipose$ID1=="PEER"] <- "PEER (111)"
master.adipose$ID1[master.adipose$ID1=="RUVCorr"] <- "RUVCorr (66)"

master.adipose$ID2[master.adipose$ID2=="None"] <- "None (79)"
master.adipose$ID2[master.adipose$ID2=="CONFETI"] <- "CONFETI (93)"
master.adipose$ID2[master.adipose$ID2=="Known\nCovariates"] <- "Known\nCovariates (102)"
master.adipose$ID2[master.adipose$ID2=="PC"] <- "PC (136)"
master.adipose$ID2[master.adipose$ID2=="PEER"] <- "PEER (111)"
master.adipose$ID2[master.adipose$ID2=="RUVCorr"] <- "RUVCorr (66)"





master.blood$ID1[master.blood$ID1=="None"] <- "None (49)"
master.blood$ID1[master.blood$ID1=="CONFETI"] <- "CONFETI (92)"
master.blood$ID1[master.blood$ID1=="Known\nCovariates"] <- "Known\nCovariates (61)"
master.blood$ID1[master.blood$ID1=="PC"] <- "PC (102)"
master.blood$ID1[master.blood$ID1=="PEER"] <- "PEER (111)"
master.blood$ID1[master.blood$ID1=="RUVCorr"] <- "RUVCorr (55)"

master.blood$ID2[master.blood$ID2=="None"] <- "None (49)"
master.blood$ID2[master.blood$ID2=="CONFETI"] <- "CONFETI (92)"
master.blood$ID2[master.blood$ID2=="Known\nCovariates"] <- "Known\nCovariates (61)"
master.blood$ID2[master.blood$ID2=="PC"] <- "PC (102)"
master.blood$ID2[master.blood$ID2=="PEER"] <- "PEER (111)"
master.blood$ID2[master.blood$ID2=="RUVCorr"] <- "RUVCorr (55)"



p1 <- ggplot(master.intestine, aes(ID1,ID2,fill=n)) + geom_tile() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_text(color="black", size=14), axis.text.y=element_text(color="black", size=14),plot.title = element_text(size=15)) + scale_fill_gradient(low='white', high="#9999FF", na.value='grey50') + geom_tile(colour="white",size=0.2)+ guides(fill=guide_legend(title="Number of module pairs\nwith Jaccard index > 0.5", title.theme=element_text(size=14)))+
labs(x="",y="",title="Small intestine")+ geom_text(aes(label = n), size=7)

p2 <- ggplot(master.spleen, aes(ID1,ID2,fill=n)) + geom_tile() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_text(color="black", size=14), axis.text.y=element_text(color="black", size=14),plot.title = element_text(size=15)) + scale_fill_gradient(low='white', high="#9999FF", na.value='grey50') + geom_tile(colour="white",size=0.2)+ guides(fill=guide_legend(title="Number of module pairs\nwith Jaccard index > 0.5", title.theme=element_text(size=14)))+
labs(x="",y="",title="Spleen")+ geom_text(aes(label = n), size=7)

p3 <- ggplot(master.heart, aes(ID1,ID2,fill=n)) + geom_tile() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_text(color="black", size=14), axis.text.y=element_text(color="black", size=14),plot.title = element_text(size=15)) + scale_fill_gradient(low='white', high="#9999FF", na.value='grey50') + geom_tile(colour="white",size=0.2)+ guides(fill=guide_legend(title="Number of module pairs\nwith Jaccard index > 0.5", title.theme=element_text(size=14)))+
labs(x="",y="",title="Heart-left ventricle")+ geom_text(aes(label = n), size=7)

p4 <- ggplot(master.adipose, aes(ID1,ID2,fill=n)) + geom_tile() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_text(color="black", size=14), axis.text.y=element_text(color="black", size=14),plot.title = element_text(size=15)) + scale_fill_gradient(low='white', high="#9999FF", na.value='grey50') + geom_tile(colour="white",size=0.2)+ guides(fill=guide_legend(title="Number of module pairs\nwith Jaccard index > 0.5", title.theme=element_text(size=14)))+
labs(x="",y="",title="Adipose")+ geom_text(aes(label = n), size=7)

p5 <- ggplot(master.muscle, aes(ID1,ID2,fill=n)) + geom_tile() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_text(color="black", size=14), axis.text.y=element_text(color="black", size=14),plot.title = element_text(size=15)) + scale_fill_gradient(low='white', high="#9999FF", na.value='grey50') + geom_tile(colour="white",size=0.2)+ guides(fill=guide_legend(title="Number of module pairs\nwith Jaccard index > 0.5", title.theme=element_text(size=14)))+
labs(x="",y="",title="Skeletal muscle")+ geom_text(aes(label = n), size=7)

p6 <- ggplot(master.blood, aes(ID1,ID2,fill=n)) + geom_tile() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_text(color="black", size=14), axis.text.y=element_text(color="black", size=14),plot.title = element_text(size=15)) + scale_fill_gradient(low='white', high="#9999FF", na.value='grey50') + geom_tile(colour="white",size=0.2)+ guides(fill=guide_legend(title="Number of module pairs\nwith Jaccard index > 0.5", title.theme=element_text(size=14)))+
labs(x="",y="",title="Whole blood") + geom_text(aes(label = n), size=7)


pdf("module_overlap_heatmaps.pdf", width=22, height=13)
pmaster <- (p1 | p2 ) / (p3 | p4) / (p5 | p6) 
pmaster
dev.off()

