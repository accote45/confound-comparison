# adapted from https://github.com/jsomekh/BCeF_

plotBCeF<-function(input.edata, input.covariates.df, input.gold.standard, input.edata.description = "", input.adjustment.method.description = "adjusted", color.to.use = c("black","red"))
{
  if(!require("pROC")){
    print("error: pROC package was not installed")
    return(0)
  }
  raw.edata           = as.matrix(input.edata)
  t.raw.edata         = t(raw.edata)

  adjusted.reads = as.matrix(input.covariates.df)
  #generate the ROC, AUC comparing raw and adjusted data
  gold.standard.to.use            = input.gold.standard
  colnames(gold.standard.to.use)  = c("Gene1","Gene2","Confidence")
  rownames(gold.standard.to.use)  = c(1:nrow(gold.standard.to.use))
  
  bins.cors.df                  = gold.standard.to.use
  #bins.cors.df$binAll.raw.r     = NA#add a zero column
  bins.cors.df$binAll.raw.pval  = NA
  #bins.cors.df$binAll.r         = NA#add a zero column
  bins.cors.df$binAll.pval      = NA
  
  for(j in 1:nrow(gold.standard.to.use))
  {
    gene1.ensmbl = gold.standard.to.use$Gene1[j]
    gene2.ensmbl = gold.standard.to.use$Gene2[j]
    
    if( (gene1.ensmbl %in% rownames(raw.edata)) & (gene2.ensmbl %in% rownames(raw.edata)))
    {
      gene1.bin.all.raw.vec   = gene2.bin.all.raw.vec = 0
      gene1.bin.all.raw.vec   = raw.edata[gene1.ensmbl,]
      gene2.bin.all.raw.vec   = raw.edata[gene2.ensmbl,]
      result.all.raw          = stats::cor.test(gene1.bin.all.raw.vec, gene2.bin.all.raw.vec, method = "pearson")
      
      gene1.bin.all.vec       = gene2.bin.all.vec = 0
      gene1.bin.all.vec       = adjusted.reads[gene1.ensmbl,]
      gene2.bin.all.vec       = adjusted.reads[gene2.ensmbl,]
      result.all              = stats::cor.test(gene1.bin.all.vec, gene2.bin.all.vec, method = "pearson")
      
      #bins.cors.df[j,"binAll.raw.r"]    = result.all.raw$estimate
      bins.cors.df[j,"binAll.raw.pval"] = result.all.raw$p.value
      #bins.cors.df[j,"binAll.r"]        = result.all$estimate 
      bins.cors.df[j,"binAll.pval"]     = result.all$p.value
    }else{
      print(paste0(j," is index of not found gene in input.edata"))
      j = j+1
    }
  }

  pvals.df                = bins.cors.df[-c(which(is.na(bins.cors.df$binAll.raw.pval))),]
  rownames(pvals.df)      = c(1:nrow(pvals.df))
  gold.standard.bool      = pvals.df[,"Confidence"]
  
  pvals.df                = cbind(pvals.df[,"binAll.raw.pval"],
                                  pvals.df[,"binAll.pval"])
  pvals.df                = apply(pvals.df, 2, function(x) p.adjust(x, method = "BH"))#multiple correction for each column
  dataset.col.names       = c("Raw", input.adjustment.method.description)
  colnames(pvals.df)      = dataset.col.names
  
  pvals.df.log            = 0
  pvals.df.log            = -log(pvals.df,10)
  tmp                     = pvals.df.log
  tmp[is.infinite(tmp)]   = 0
  max.val                 = max(tmp)
  pvals.df.log[is.infinite(pvals.df.log)]<-(max.val+10)#change Inf to max value
  #pvals.df.log            = pvals.df.log[, c(input.adjustment.method.description, "Raw")]
  legend.to.use           = dataset.col.names
  
  cases.gold              = length(which(gold.standard.bool==1))
  control.gold            = length(which(gold.standard.bool==0))
  roc_obj_linear          = pROC::roc(gold.standard.bool, pvals.df.log[,"Raw"])
  l.auc                   = pROC::auc(roc_obj_linear)
  auc.vec                 <-l.auc
  plot(roc_obj_linear, 
       main     = paste0(input.edata.description, " (gold standard cases/controls are ",cases.gold,"/",control.gold,")"),
       col      = color.to.use[1],
       cex.main = 0.8)
  
  final.legend.to.use = paste0(legend.to.use[1], " (auc=",round(l.auc,3),")")
  
  for(i in 2:ncol(pvals.df.log))
  {
    roc_obj_tmp   = pROC::roc(gold.standard.bool, pvals.df.log[,i])
    tmp.auc       = pROC::auc(roc_obj_tmp)
    lines(roc_obj_tmp, col = color.to.use[i])
    final.legend.to.use    = c(final.legend.to.use, paste0(legend.to.use[i], " (auc=",round(tmp.auc,3),")"))
  }
  
  legend("bottomright",
         cex    = 1,
         legend = final.legend.to.use, 
         col    = color.to.use, 
         lwd    = 0.8, 
         lty    = c(1,1,1,1))#type of line
}



pdf("tissue_network.pdf")
plotBCeF(input.edata=dat,input.covariates.df=adj.dat,input.gold.standard=new)
dev.off()
