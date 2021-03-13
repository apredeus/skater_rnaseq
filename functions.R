pca12t = function (exp_table,cond_table,color,labels,title) {
  library(ggplot2)
  library(ggrepel)
  
  ## counts should be log2-transformed and normalized
  pca         <- prcomp(t(exp_table))
  var_PC1     <- round(summary(pca)$importance[2,1],3)
  var_PC2     <- round(summary(pca)$importance[2,2],3)
  sumvar      <- var_PC1+var_PC2
  
  pcat        <- pca$x[,1:2]
  pcat_ann    <- merge(pcat,cond_table,by=0)
  
  xlab  <- paste("PC1, ",var_PC1*100,"% variance explained",sep="")
  ylab  <- paste("PC2, ",var_PC2*100,"% variance explained",sep="")
  
  ggplot(pcat_ann, aes(PC1, PC2)) + geom_point(aes(color = pcat_ann[,colnames(pcat_ann) %in% color]),size=3) + 
    geom_text_repel(data=pcat_ann,aes(label=pcat_ann[,colnames(pcat_ann) %in% labels])) +
    ggtitle(title) + labs(x=xlab,y=ylab) + scale_colour_discrete(name=color) + theme_bw()
}

### now we make a function to automatically subset, normalize, assess, and do DE on these datasets 
limma_donor_correct = function(gse,cond,gpl_ann,log2transformed = T,normalized = T) {
  library(limma)
  library(sva)
  library(ggplot2)
  library(patchwork)
  
  ## select only needed samples; probes that are genes; log2 transform, if needed 
  exp      <- exprs(gse)
  exp.filt <- exp[,colnames(exp) %in% rownames(cond)]
  
  ## log2-transform if it wasn't 
  if (! log2transformed) {
    cat(sprintf("Matrix will be log2-transformed!\n"))
    exp.filt <- log2(exp.filt + 1)
  }
  
  ## normalize 
  if (! normalized) {
    cat(sprintf("Matrix will be normalized using limma's normalizeBetweenArrays!\n"))
    exp.filt <- normalizeBetweenArrays(exp.filt,method = "quantile")
  }
  
  exp.ann  <- merge(gpl_ann[,c(2,4)],exp.filt,by = 0)
  exp.af   <- exp.ann[exp.ann$Gene_name != "NONE",]
  rownames(exp.af) <- exp.af$Row.names
  exp.af$Row.names <- NULL
  
  ## take top 6k genes, do ComBat correction, compare PCA plots
  ann              <- exp.af[,c(1,2)]
  exp.sorted       <- exp.af[order(rowMeans(exp.af[,3:ncol(exp.af)]), decreasing = T),]
  exp.6k           <- head(exp.sorted, 6000)[,3:ncol(exp.af)]
  mod.6k           <- model.matrix(~ exercise, cond)
  cbat.6k          <- ComBat(exp.6k, batch = cond$donor, mod = mod.6k)
  
  ## make PCA plots
  p1 <- pca12t(exp.6k,cond,"donor","Row.names","Uncorrected PCA (col. by donor)")
  p2 <- pca12t(exp.6k,cond,"exercise","Row.names","Uncorrected PCA (col. by exercise)")
  p3 <- pca12t(cbat.6k,cond,"donor","Row.names","ComBat corrected PCA (col. by donor)") 
  p4 <- pca12t(cbat.6k,cond,"exercise","Row.names","ComBat corrected PCA (col. by exercise)") 

  ## limma differential expression with donor effect
  Dn <- cond$donor
  Ex <- factor(cond$exercise,levels=c("before","after"))
  design <- model.matrix(~ Dn + Ex)
  fit <- lmFit(exp.af[,3:ncol(exp.af)], design)
  fit <- eBayes(fit)
  de <- topTable(fit, n=100000, coef="Exafter")
  de.ann <- merge(ann,de,by = 0)
  
  ## plot t-values boxplot
  p5 <- ggplot(de.ann,aes(y=t)) + geom_boxplot() + theme_bw() + ggtitle("t statistic distribution")
  de.ann$Sign <- ifelse(de.ann$adj.P.Val <= 0.1 & de.ann$logFC > 0, "Up", "Not sig")
  de.ann$Sign <- ifelse(de.ann$adj.P.Val <= 0.1 & de.ann$logFC < 0, "Dn", de.ann$Sign)
  p6 <- ggplot(de.ann,aes(x=AveExpr,y=logFC,color=Sign)) + geom_point(size=0.4) + 
    scale_color_manual(values=c("Not sig" = "#999999", "Up" = "#FF0000", "Dn" = "#0000FF")) + theme_bw() + ggtitle("MA plot, FDR 10% sign cutoff")
  
  print((p1 + p2)/(p3 + p4)/(p5 + p6))
  return(de.ann)
}

### use DE output from limma to run fgsea

run_fgsea = function(de,pathways,protein_coding,ds.name) {
  library(fgsea)
  library(tidyr)
  library(BiocParallel)
  
  ## select Gene_name and t-stat columns from the limma output
  de.rnk         <- de[,c(2,6)]
  ## split rows that have multiple genes separated by comma
  de.split       <- as.data.frame(tidyr::separate_rows(de.rnk,"Gene_name",sep = ","))
  ## only select protein coding genes
  de.filt        <- de.split[de.split$Gene_name %in% protein_coding,]
  ## keep one row with highest abs. t stat per gene name 
  de.agg         <- aggregate(t ~ Gene_name, data = de.filt, function(x) ifelse(mean(x)>0,max(x),min(x)))
  ## prep the vector for fGSEA reqs 
  cat(sprintf("Row counts in dataset %s: initial diffexp: %d, split: %d, protein only: %d, collapsed: %d",
              ds.name,dim(de.rnk)[1],dim(de.split)[1],dim(de.filt)[1],dim(de.agg)[1]))
  rnk            <- setNames(de.agg$t,toupper(de.agg$Gene_name))
  fgsea.out      <- fgsea(pathways,rnk,minSize = 15,maxSize = 500,eps = 0.0,nproc = 4)
  fgsea.out      <- as.data.frame(fgsea.out[,c(1,6)])
  colnames(fgsea.out)[2] <- ds.name
  return(fgsea.out)
}
