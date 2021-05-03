

files <- list.files(pattern = "unsigned_pearson_min_mod_size10connectivity.txt$", recursive = TRUE, full.names=TRUE)
files <- files[!grepl("skin",files)]
files <- files[!grepl("artery",files)]
files <- files[!(files %like% "null")]
listOfFiles <- lapply(files, function(x) read.table(x, header = TRUE)) 


files2 <- list.files(pattern = "unsigned_pearson_min_mod_size10module_assigned_master.txt$", recursive = TRUE, full.names=TRUE)
files2 <- files2[!grepl("skin",files2)]
files2 <- files2[!grepl("artery",files2)]
files2 <- files2[!(files2 %like% "null")]
listOfFiles2 <- lapply(files2, function(x) read.table(x, header = TRUE)) 

# add module assignment to connectivity master file
listOfFiles <- lapply(seq_along(listOfFiles), function(i) cbind(listOfFiles[[i]],gene.names=rownames(listOfFiles[[i]])))

connect.dat <- lapply(seq_along(listOfFiles), function(i) merge(listOfFiles[[i]],listOfFiles2[[i]], by="gene.names"))
connect.dat <- lapply(connect.dat,function(x) x[x$mergedColors!="grey",])

names(connect.dat) <- gsub("unsigned_pearson_min_mod_size10numbergenes_permodule.txt|/|\\.","",files)
 
# calculate network density per module

master <- lapply(seq_along(connect.dat),function(i) aggregate(connect.dat[[i]][,3],list(connect.dat[[i]]$mergedColors),sum))
master <- lapply(master,setNames,c("mergedColors","x"))

files3 <- list.files(pattern = "numbergenes_permodule.txt$", recursive = TRUE, full.names=TRUE)
files3 <- files3[!grepl("skin",files3)]
files3 <- files3[!grepl("artery",files3)]
files3 <- files3[!(files3 %like% "null")]


listOfFiles3 <- lapply(files3, function(x) read.table(x, header = TRUE)) 
# remove grey module
listOfFiles3 <- lapply(listOfFiles3,function(x) x[x$mergedColors!="grey",])

master.fin <- lapply(seq_along(listOfFiles3), function(i) merge(listOfFiles3[[i]],master[[i]], by="mergedColors"))

for (i in seq_along(master.fin)){
  master.fin[[i]]$density <- master.fin[[i]]$x/(master.fin[[i]]$N*master.fin[[i]]$N)
}

names(master.fin) <- gsub("unsigned_pearson_min_mod_size10connectivity.txt|/|\\.","",files)


master.fin.dat <- lapply(seq_along(master.fin), function(i) cbind(master.fin[[i]],type=names(master.fin[i])))
names(master.fin.dat) <- gsub("unsigned_pearson_min_mod_size10connectivity.txt|/|\\.","",files)

master.fin.dat <- ldply(master.fin.dat,rbind)
master.fin.dat$type <- gsub("small_intestine","intestine",master.fin.dat$type)
master.fin.dat <- separate(master.fin.dat, type, into=c("Tissue","Adjustment"),sep="_")
master.fin.dat$mergedColors <- NULL

master.fin.dat$Adjustment[master.fin.dat$Adjustment=="known"] <- "Known\nCovariates"
master.fin.dat$Adjustment[master.fin.dat$Adjustment=="basic"] <- "None"
master.fin.dat$Adjustment[master.fin.dat$Adjustment=="ruvcorr"] <- "RUVCorr"
master.fin.dat$Adjustment[master.fin.dat$Adjustment=="peer"] <- "PEER"
master.fin.dat$Adjustment[master.fin.dat$Adjustment=="pc"] <- "PC"
master.fin.dat$Adjustment[master.fin.dat$Adjustment=="confeti"] <- "CONFETI"

master.fin.dat$Tissue[master.fin.dat$Tissue=="adipose"] <- "Adipose"
master.fin.dat$Tissue[master.fin.dat$Tissue=="blood"] <- "Whole blood"
master.fin.dat$Tissue[master.fin.dat$Tissue=="heart"] <- "Heart-left ventricle"
master.fin.dat$Tissue[master.fin.dat$Tissue=="muscle"] <- "Skeletal muscle"
master.fin.dat$Tissue[master.fin.dat$Tissue=="intestine"] <- "Small intestine"
master.fin.dat$Tissue[master.fin.dat$Tissue=="spleen"] <- "Spleen"

master.fin.dat$Adjustment <- factor(master.fin.dat$Adjustment, levels=c("None","Known\nCovariates","RUVCorr","PC","PEER","CONFETI"))

pdf("mod_properties_density.pdf", width=9, height=8)
p <- ggplot(master.fin.dat, aes(x=Adjustment, y=density, fill=Tissue)) + geom_boxplot(outlier.shape=NA) + theme_bw() + ylab("Module Density") + ylim(limits = c(0, .25)) +
	scale_fill_manual(values=jjsparkles_palette) + theme(axis.text.x = element_text(color = "black", size = 13, angle = 0, , vjust=0.5, face = "plain"),
		 axis.text.y = element_text(color = "black", size = 13, angle = 0, , vjust=0.5, face = "plain"),
                axis.title.y=element_text(color = "black", size = 17, angle = 90, , vjust=0.5, face = "plain"),
                 axis.title.x=element_text(color = "black", size = 17, angle = 0, , vjust=0.5, face = "plain"),
                 legend.text=element_text(size=13), legend.title=element_text(size=15))
p
dev.off()
