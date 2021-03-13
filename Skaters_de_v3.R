### analysis v3. Differential expression of RNA-seq, pathway 
## enrichment (fGSEA), and comparison with microarray experiments
## put rsem tsv files into the dir you're working in. 


library(rapportools)
library(data.table)
library(DESeq2)
library(limma)
library(fgsea)
library(BiocParallel)
library(ggplot2)
library(ggrepel)
library(reshape)
library(sva)
library(vsn)
library(DT)
library(GEOquery)
library(dplyr)
library(clusterProfiler)
library(enrichplot)

setwd("~/tmp")
source("~/tmp/functions.R") 

################################## RNA-seq differential expression analysis #########################################

## PAR genes were removed from RSEM reference 
rsem_exp      <- read.table("rsem.genes.counts.tsv", header=T, row.names = 1, check.names=F)
rsem_tpm      <- read.table("rsem.genes.TPM.tsv", header=T, row.names = 1, check.names=F)

dim(rsem_exp)
head(rsem_exp)
head(rsem_tpm)

cond          <- read.table("Conditions.txt", row.names = 1, header = T, check.names = F,stringsAsFactors = T)
ann           <- rsem_exp[,c(1:2)]
exp           <- rsem_exp[,c(3:16)]
colSums(exp[, c(1:14)])

## make 12k table for visualization, do PCA on the rlog-transformed matrix, before and after donor correction with ComBat

exp_sorted      <- exp[order(rowMeans(exp), decreasing = T),]
exp_top12k      <- head(exp_sorted, 12000)
dds_top12k      <- DESeqDataSetFromMatrix(countData=round(exp_top12k,0),colData=cond,design = ~ Donor + Cond)
mod12k          <- model.matrix(~ Cond, cond)
rlog_top12k     <- assay(rlog(dds_top12k,blind = F,fitType = "mean")) 
cbat_top12k     <- ComBat(rlog_top12k, batch = cond$Donor, mod = mod12k)
boxplot(rlog_top12k)
boxplot(cbat_top12k)

pca12t(rlog_top12k,cond,"Cond","Donor","PCA before donor effect correction")
pca12t(cbat_top12k,cond,"Cond","Donor","PCA after effect correction") 

## Let's save annotated, ComBat-corrected table for heatmap visualizations

top12k_ann      <- merge(ann, cbat_top12k, by = 0)
colnames(top12k_ann)[1] <- c("Ensembl")
write.table(top12k_ann,"Skaters_top12k_cbat.tsv",quote = F,sep = "\t",row.names = F)

############# now let's do DESeq2 analysis and make some plots for diagnostics and illustration

dds_full      <- DESeqDataSetFromMatrix(countData=round(exp,0),colData = cond,design = ~ Donor + Cond)
dds_full$Cond <- relevel(dds_full$Cond, ref = "before") ## all our log2FC are = (after/before)
dds_full      <- DESeq(dds_full)
res_full      <- results(dds_full, contrast=c("Cond","after","before"))
res_df        <- as.data.frame(res_full)
res_ann       <- merge(ann,res_df,by = 0)
res_ann$lfcSE <- NULL
colnames(res_ann)[1] <- c("Ensembl")
write.table(res_ann,"Skaters_full_diffexp.tsv",quote = F,sep = "\t",row.names = F)

## make MA plot
res_fann <- res_full[complete.cases(res_full),]
DESeq2::plotMA(res_full, main = "MA plot, before vs. after exercise",ylim = c(-4, 4),
               colLine = "black", colSig = "red",colNonSig = "gray80")

## annotate diffexp results with TPM, make filtered sets

rsem_tpm$Mean <- rowMeans(rsem_tpm[,3:16])
mean_tpm <- rsem_tpm[,17,drop = F]
colnames(mean_tpm) <- "Mean_TPM"
res_tpm <- merge(res_fann,mean_tpm,by.x = "Ensembl",by.y = "row.names")

## make Volcano plot with dot size proportional to expression
res_label <- subset(res_tpm, abs(log2FoldChange) > 2 | padj < 1e-8)

ggplot(res_tpm, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant, size = log10(Mean_TPM+1)), stroke = 0) + 
  geom_text_repel(data = res_label, aes(label = Symbol), size = 2) + 
  scale_size_continuous(range = c(0.2,3)) + 
  scale_x_continuous(limits = c(-3.5, 3.5)) + 
  scale_color_manual(values = c("red", "gray80")) + 
  xlab("log fold change") + 
  ylab("-log10 (adjusted p-value)") + 
  ggtitle("Volcano plot, before vs. after exercise") + theme_bw() + 
  theme(panel.grid.major = element_line(size = 0.25))

## make bar chart of differentially expressed gene numbers, filtered by mean TPM
cuts <- c(10,100,1000)
up <- res_tpm[res_tpm$log2FoldChange > 0 & res_tpm$padj <= 0.1,]
dn <- res_tpm[res_tpm$log2FoldChange < 0 & res_tpm$padj <= 0.1,]
de_genes <- data.frame(nrow(up),nrow(dn))
colnames(de_genes) <- c("up","dn")

for (i in cuts) {
  up <- res_tpm[res_tpm$log2FoldChange > 0 & res_tpm$padj <= 0.1 & res_tpm$Mean_TPM > i,]
  dn <- res_tpm[res_tpm$log2FoldChange < 0 & res_tpm$padj <= 0.1 & res_tpm$Mean_TPM > i,]
  cat(sprintf("Cutoff: %d, up: %d, down: %d\n",i,nrow(up),nrow(dn)))
  de_genes <- rbind(de_genes,c(nrow(up),nrow(dn)))
}

de_genes$cutoff <- c("TPM > 0", "TPM > 10", "TPM > 100", "TPM > 1000")
de_melt <- reshape2::melt(de_genes)

ggplot(de_melt, aes(x=cutoff, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.5)  + 
  scale_fill_manual(values = c("red", "blue")) + 
  xlab("Mean expression cutoff (TPM)") + 
  ylab("Differentially expressed genes") + 
  ggtitle("Differentially expressed genes, before vs. after exercise") + theme_bw() + 
  theme(panel.grid.major = element_line(size = 0.25))

## Save filtered differentially expressed gene sets

write.table(res_tpm[res_tpm$padj <= 0.1 & res_tpm$Gene_type == "protein_coding",],
            "Skaters_all_de.tsv",quote = F,sep = "\t",row.names = F)
write.table(res_tpm[res_tpm$padj <= 0.1 & res_tpm$Mean_TPM >= 10 & res_tpm$Gene_type == "protein_coding",],
            "Skaters_all_de_tpm10.tsv",quote = F,sep = "\t",row.names = F)
write.table(res_tpm[res_tpm$padj <= 0.1 & res_tpm$Mean_TPM >= 100 & res_tpm$Gene_type == "protein_coding",],
            "Skaters_all_de_tpm100.tsv",quote = F,sep = "\t",row.names = F)

################################## Pathway enrichment with enrichR and fGSEA #######################################

## make lists of protein-coding up- and down-regulated genes (as gene symbols)
up_symb <- res_tpm[res_tpm$log2FoldChange > 0 & res_tpm$padj <= 0.1 & res_tpm$Gene_type == "protein_coding",]$Symbol
dn_symb <- res_tpm[res_tpm$log2FoldChange < 0 & res_tpm$padj <= 0.1 & res_tpm$Gene_type == "protein_coding",]$Symbol
length(up_symb)
length(dn_symb)

### GMT files can be downloaded from http://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb
h.all          <- gmtPathways("~/tmp/h.all.v7.2.symbols.gmt")
c2.cp          <- gmtPathways("~/tmp/c2.cp.v7.2.symbols.gmt")
c7.all         <- gmtPathways("~/tmp/c7.all.v7.2.symbols.gmt")

## make term-2-gene data frames for clusterProfiler
h_t2g <- reshape2::melt(h.all)
colnames(h_t2g) <- c("Symbol","Pathway")
h_t2g <- h_t2g[,c(2,1)]

cp_t2g <- reshape2::melt(c2.cp)
colnames(cp_t2g) <- c("Symbol","Pathway")
cp_t2g <- cp_t2g[,c(2,1)]

c7_t2g <- reshape2::melt(c7.all)
colnames(c7_t2g) <- c("Symbol","Pathway")
c7_t2g <- c7_t2g[,c(2,1)]

## calculate and plot enrichments using clusterProfiler 
universe <- res_ann[res_ann$Gene_type == "protein_coding",]$Symbol
enr_up_h <- enricher(gene = up_symb, TERM2GENE = h_t2g)
enr_dn_h <- enricher(gene = dn_symb, TERM2GENE = h_t2g, universe = res_ann$Symbol)
enr_up_cp <- enricher(gene = up_symb, TERM2GENE = cp_t2g, universe = universe)
enr_dn_cp <- enricher(gene = dn_symb, TERM2GENE = cp_t2g, universe = universe)

dotplot(enr_up_h, showCategory = 10, title = "Gene overlap, up-regulated genes vs. MsigDB H")
dotplot(enr_up_cp, showCategory = 10,title = "Gene overlap, up-regulated genes vs. MsigDB C2:CP")
dotplot(enr_dn_cp, showCategory = 10,title = "Gene overlap, down-regulated genes vs. MsigDB C2:CP")

## run similar enrichments using fGSEA
rnk            <- aggregate(stat ~ Symbol,data = res_tpm[res_tpm$Gene_type == "protein_coding",], function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk            <- setNames(rnk$stat,toupper(rnk$Symbol))
length(rnk)

gsea.h         <- fgsea(h.all,rnk,minSize = 10,maxSize = 500,eps = 0.0,nproc = 4)
gsea.cp        <- fgsea(c2.cp,rnk,minSize = 10,maxSize = 500,eps = 0.0,nproc = 4)

gsea_up_h <- gsea.h[gsea.h$NES > 0 & gsea.h$padj < 0.05,]
gsea_up_h <- gsea_up_h[order(gsea_up_h$padj),]
gsea_up_h <- gsea_up_h[1:10,]
gsea_up_h$Overlap <- lengths(gsea_up_h$leadingEdge)

gsea_dn_h <- gsea.h[gsea.h$NES < 0 & gsea.h$padj < 0.05,] 
gsea_dn_h <- gsea_up_h[order(gsea_up_h$padj),] ## very similar situation, only MYC v2 is very significant
gsea_up_cp <- gsea.cp[gsea.cp$NES > 0 & gsea.cp$padj < 0.05,]
gsea_up_cp <- gsea_up_cp[order(gsea_up_cp$padj),]
gsea_up_cp <- gsea_up_cp[1:10,]
gsea_up_cp$Overlap <- lengths(gsea_up_cp$leadingEdge)

gsea_dn_cp <- gsea.cp[gsea.cp$NES < 0 & gsea.cp$padj < 0.05,] 
gsea_dn_cp <- gsea_dn_cp[order(gsea_dn_cp$padj),] 
gsea_dn_cp <- gsea_dn_cp[1:10,]
gsea_dn_cp$Overlap <- lengths(gsea_dn_cp$leadingEdge)

## let's plot same enrichments as above but for GSEA
ggplot(gsea_up_h,aes(x = NES,y = factor(pathway, levels = gsea_up_h$pathway[order(gsea_up_h$NES)]),size = Overlap,color = -log(padj))) + 
  geom_point() + scale_color_gradient(low = "blue",high = "red") + 
  xlab("Normalizes enrichment score (NES)") + 
  ylab("MsigDB hallmark (H) pathways") + 
  ggtitle("GSEA enrichment, up-regulated pathways from MsigDB H") + theme_bw()

ggplot(gsea_up_cp,aes(x = NES,y = factor(pathway, levels = gsea_up_cp$pathway[order(gsea_up_cp$NES)]),size = Overlap,color = -log(padj))) + 
  geom_point() + scale_color_gradient(low = "blue",high = "red") + 
  xlab("Normalizes enrichment score (NES)") + 
  ylab("MsigDB canonical (CP) pathways") + 
  ggtitle("GSEA enrichment, up-regulated pathways from MsigDB CP") + theme_bw()

ggplot(gsea_dn_cp,aes(x = -NES,y = factor(pathway, levels = gsea_dn_cp$pathway[order(-gsea_dn_cp$NES)]),size = Overlap,color = -log(padj))) + 
  geom_point() + scale_color_gradient(low = "blue",high = "red") + 
  xlab("Negative normalizes enrichment score (-NES)") + 
  ylab("MsigDB canonical (CP) pathways") + 
  ggtitle("GSEA enrichment, down-regulated pathways from MsigDB CP") + theme_bw()


################################## Put our dataset in the context of other datasets ########################


## use GEOquery package to get the selected datasets
m2014.gse     <- getGEO("GSE51216",GSEMatrix = T)[[1]] ## GPL6480, log2-tr, normalized
t2013.gse     <- getGEO("GSE46075",GSEMatrix = T)[[1]] ## GPL6244, not log2, normalized
b2007.gse     <- getGEO("GSE3606",GSEMatrix = T)[[1]] ## GPL571, not log2, normalized
c2004.gse     <- getGEO("GSE1140",GSEMatrix = T)[[1]] ## GPL96, not log2, normalized
s2012.gse     <- getGEO("GSE28498",GSEMatrix = T)[[1]] ## GPL6244, not log2, normalized (t2 & t3 have clear batch - probably done separately)
ra2009a.gse   <- getGEO("GSE14642",GSEMatrix = T)[[1]] ## GPL570, not log2, normalized
ra2009b.gse   <- getGEO("GSE11761",GSEMatrix = T)[[1]] ## GPL570, not log2, normalized
n2010.gse     <- getGEO("GSE18966",GSEMatrix = T)[[1]] ## GPL4133, not log2, kind of normalized (?)

## custom selection of samples, and defining conditions
### Mukerjee 2014, GSE 51216
gsms <- c("GSM1240621","GSM1240622","GSM1240623","GSM1240624","GSM1240625","GSM1240626","GSM1240627","GSM1240628","GSM1240629")
exercise <- c("before","before","before","before","before","after","after","after","after")
donor <- c("d1","d2","d2","d3","d4","d1","d2","d3","d4")
m2014.cond <- as.data.frame(cbind(exercise,donor),stringsAsFactors = T)
rownames(m2014.cond) <- gsms

## Tonevitsky 2013, GSE46075 - only T1 and T2 
gsms <- c("GSM1122850","GSM1122857","GSM1122864","GSM1122878","GSM1122885","GSM1128806","GSM1128813","GSM1128820",
          "GSM1122851","GSM1122858","GSM1122865","GSM1122879","GSM1122886","GSM1128807","GSM1128814","GSM1128821")
exercise <- c("before","before","before","before","before","before","before","before",
              "after","after","after","after","after","after","after","after")
donor <- c("d1","d2","d3","d4","d5","d6","d7","d8",
           "d1","d2","d3","d4","d5","d6","d7","d8")
t2013.cond <- as.data.frame(cbind(exercise,donor),stringsAsFactors = T)
rownames(t2013.cond) <- gsms

## Buttner 2007, GSE3606 - only exhaustive exercise 
gsms <- c("GSM82670","GSM82671","GSM82674","GSM82675","GSM82676","GSM82680","GSM82682","GSM82685","GSM82686","GSM82689")
exercise <- c("before","after","before","before","after","after","before","before","after","after")
donor <- c("d1","d1","d2","d3","d3","d2","d4","d5","d5","d4")
b2007.cond <- as.data.frame(cbind(exercise,donor),stringsAsFactors = T)
rownames(b2007.cond) <- gsms

## Connely 2004, GSE1140 - only before and after (A & B, but not C)
gsms <- c("GSM19083","GSM19084","GSM19085","GSM19087","GSM19089","GSM19090","GSM19092","GSM19093","GSM19095","GSM19096")
exercise <- c("before","before","after","after","before","after","before","after","before","after")
donor <- c("d1","d2","d1","d2","d3","d3","d4","d4","d5","d5")
c2004.cond <- as.data.frame(cbind(exercise,donor),stringsAsFactors = T)
rownames(c2004.cond) <- gsms

## Sakharov 2012, GSE28498 - only t0 and t1
gsms <- c("GSM706084","GSM706085","GSM706086","GSM706087","GSM706088","GSM706089","GSM706090","GSM706091","GSM706092",
          "GSM706093","GSM706094","GSM706095","GSM706096","GSM706097","GSM706098","GSM706099","GSM706100","GSM706101",
          "GSM706102","GSM706103","GSM706104","GSM706105","GSM706106","GSM706107","GSM706108","GSM706109","GSM706110",
          "GSM706111","GSM706112","GSM706113","GSM706114","GSM706115","GSM706116","GSM706117","GSM706118","GSM706119","GSM706120","GSM706121")
exercise <- c("before","after","before","after","before","after","before","after","before","after","before","after",
              "before","after","before","after","before","after","before","after","before","after","before","after",
              "before","after","before","after","before","after","before","after","before","after","before","after","before","after")
donor <- c("d1","d1","d2","d2","d3","d3","d4","d4","d5","d5","d6","d6","d7","d7","d8","d8","d9","d9","d10","d10",
           "d11","d11","d12","d12","d13","d13","d14","d14","d15","d15","d16","d16","d17","d17","d18","d18","d19","d19")
s2012.cond <- as.data.frame(cbind(exercise,donor),stringsAsFactors = T)
rownames(s2012.cond) <- gsms

## Radom-Azik 2009a, GSE14642 (girls), all experiments
gsms <- paste("GSM365",as.character(496:535),sep = "")
exercise <- rep(c("before","after"),times = 20)
donor <- rep(paste("d",as.character(1:20),sep = ""),each = 2)
ra2009a.cond <- as.data.frame(cbind(exercise,donor),stringsAsFactors = T)
rownames(ra2009a.cond) <- gsms

## Radom-Azik 2009b GSE11761 (boys), all experiments
gsms <- paste("GSM2978",as.character(52:91),sep = "")
exercise <- rep(c("before","after"),times = 20)
donor <- rep(paste("d",as.character(1:20),sep = ""),each = 2)
ra2009b.cond <- as.data.frame(cbind(exercise,donor),stringsAsFactors = T)
rownames(ra2009b.cond) <- gsms

## Nakamura 2010, GSE18966, before and after 4h exercise (0, 4h)
gsms <- c("GSM469524","GSM469525","GSM469528","GSM469529","GSM469532","GSM469533","GSM469536","GSM469537","GSM469540","GSM469541")
exercise <- c("before","after","before","after","before","after","before","after","before","after")
donor <- rep(paste("d",as.character(1:5),sep = ""),each = 2)
n2010.cond <- as.data.frame(cbind(exercise,donor),stringsAsFactors = T)
rownames(n2010.cond) <- gsms

## custom annotations are used with updated gene names etc
gpl96.ann    <- read.table("GPL96.ann.tsv",header = T, row.names = 1)
gpl570.ann   <- read.table("GPL570.ann.tsv",header = T, row.names = 1)
gpl571.ann   <- read.table("GPL571.ann.tsv",header = T, row.names = 1)
gpl6244.ann  <- read.table("GPL6244.ann.tsv",header = T, row.names = 1)
gpl6480.ann  <- read.table("GPL6480.ann.tsv",header = T, row.names = 1)
gpl4133.ann  <- read.table("GPL4133.ann.tsv",header = T, row.names = 1)

## 
m2014.de     <- limma_donor_correct(m2014.gse,m2014.cond,gpl6480.ann,log2transformed = T,normalized = T)
t2013.de     <- limma_donor_correct(t2013.gse,t2013.cond,gpl6244.ann,log2transformed = F,normalized = F) ## groups look weird
b2007.de     <- limma_donor_correct(b2007.gse,b2007.cond,gpl571.ann,log2transformed = F,normalized = T) 
c2004.de     <- limma_donor_correct(c2004.gse,c2004.cond,gpl96.ann,log2transformed = F,normalized = T)
s2012.de     <- limma_donor_correct(s2012.gse,s2012.cond,gpl6244.ann,log2transformed = F,normalized = T)
ra2009a.de   <- limma_donor_correct(ra2009a.gse,ra2009a.cond,gpl570.ann,log2transformed = F,normalized = T)
ra2009b.de   <- limma_donor_correct(ra2009b.gse,ra2009b.cond,gpl570.ann,log2transformed = F,normalized = T)
n2010.de     <- limma_donor_correct(n2010.gse,n2010.cond,gpl4133.ann,log2transformed = F,normalized = F)

##  most of these look amazingly consistent. 

## let's do fGSEA vs hallmark and then PCA 

h.all          <- gmtPathways("~/tmp/h.all.v7.2.symbols.gmt")
protein_coding <- ann[ann$Gene_type == "protein_coding",]$Symbol


m2014.fgsea    <- run_fgsea(m2014.de,h.all,protein_coding,"m2014")
b2007.fgsea    <- run_fgsea(b2007.de,h.all,protein_coding,"b2007")
c2004.fgsea    <- run_fgsea(c2004.de,h.all,protein_coding,"c2004")
s2012.fgsea    <- run_fgsea(s2012.de,h.all,protein_coding,"s2012")
ra2009a.fgsea  <- run_fgsea(ra2009a.de,h.all,protein_coding,"ra2009a")
ra2009b.fgsea  <- run_fgsea(ra2009b.de,h.all,protein_coding,"ra2009b")
n2010.fgsea    <- run_fgsea(n2010.de,h.all,protein_coding,"n2010")

our.h <- as.data.frame(gsea.h[,c(1,6)])
colnames(our.h)[2] <- "ski"

## extremely useful function
all_fgsea = Reduce(function(...) merge(..., by = "pathway"), 
                   list(our.h,m2014.fgsea,b2007.fgsea,c2004.fgsea,s2012.fgsea,ra2009a.fgsea,ra2009b.fgsea,n2010.fgsea))
rownames(all_fgsea) <- all_fgsea$pathway
all_fgsea$pathway <- NULL

pca         <- prcomp(t(all_fgsea))
var_PC1     <- round(summary(pca)$importance[2,1],3)
var_PC2     <- round(summary(pca)$importance[2,2],3)
sumvar      <- var_PC1+var_PC2
pcat        <- as.data.frame(pca$x[,1:2])

xlab  <- paste("PC1, ",var_PC1*100,"% variance explained",sep="")
ylab  <- paste("PC2, ",var_PC2*100,"% variance explained",sep="")

ggplot(pcat, aes(PC1, PC2)) + geom_point(size=3) + geom_text_repel(label=rownames(pcat)) + 
  ggtitle("Hallmark fgsea PCA") + labs(x=xlab,y=ylab) + theme_bw()

##### try the same for cp.

c2.cp          <- gmtPathways("~/tmp/c2.cp.v7.2.symbols.gmt")
protein_coding <- ann[ann$Gene_type == "protein_coding",]$Symbol


m2014.fgsea    <- run_fgsea(m2014.de,c2.cp,protein_coding,"m2014")
b2007.fgsea    <- run_fgsea(b2007.de,c2.cp,protein_coding,"b2007")
c2004.fgsea    <- run_fgsea(c2004.de,c2.cp,protein_coding,"c2004")
s2012.fgsea    <- run_fgsea(s2012.de,c2.cp,protein_coding,"s2012")
ra2009a.fgsea  <- run_fgsea(ra2009a.de,c2.cp,protein_coding,"ra2009a")
ra2009b.fgsea  <- run_fgsea(ra2009b.de,c2.cp,protein_coding,"ra2009b")
n2010.fgsea    <- run_fgsea(n2010.de,c2.cp,protein_coding,"n2010")

our.cp <- as.data.frame(gsea.cp[,c(1,6)])
colnames(our.cp)[2] <- "ski"

## extremely useful function
all_fgsea2 = Reduce(function(...) merge(..., by = "pathway"), 
                    list(our.cp,m2014.fgsea,b2007.fgsea,c2004.fgsea,s2012.fgsea,ra2009a.fgsea,ra2009b.fgsea,n2010.fgsea))
rownames(all_fgsea2) <- all_fgsea2$pathway
all_fgsea2$pathway <- NULL

pca         <- prcomp(t(all_fgsea2))
var_PC1     <- round(summary(pca)$importance[2,1],3)
var_PC2     <- round(summary(pca)$importance[2,2],3)
sumvar      <- var_PC1+var_PC2
pcat        <- as.data.frame(pca$x[,1:2])

xlab  <- paste("PC1, ",var_PC1*100,"% variance explained",sep="")
ylab  <- paste("PC2, ",var_PC2*100,"% variance explained",sep="")

ggplot(pcat, aes(PC1, PC2)) + geom_point(size=3) + geom_text_repel(label=rownames(pcat)) + 
  ggtitle("Hallmark fgsea PCA") + labs(x=xlab,y=ylab) + theme_bw()

## so probably should do something different 
our.de <- res_ann[,c(2,6)]
our.de <- our.de[complete.cases(our.de),]
colnames(our.de) <- c("Gene_name","ski")
all <- merge(our.de,m2014.de,by = "Gene_name")[,c(1,2,7)]
colnames(all)[3] <- "m2014"
all2 <- merge(all,b2007.de,by = "Gene_name")[,c(1:3,8)]
colnames(all2)[4] <- "b2007"
all <- merge(all,c2004.de,by = "Gene_name")[,c(1:4,9)]
colnames(all)[5] <- "c2004"
all <- merge(all,s2012.de,by = "Gene_name")[,c(1:5,10)]
colnames(all)[6] <- "s2012"
all <- merge(all,ra2009a.de,by = "Gene_name")[,c(1:6,11)]
colnames(all)[7] <- "ra2009a"
all <- merge(all,ra2009b.de,by = "Gene_name")[,c(1:7,12)]
colnames(all)[8] <- "ra2009b"
all <- merge(all,n2010.de,by = "Gene_name")[,c(1:8,13)]
colnames(all)[9] <- "n2010"



m2014.fgsea    <- run_fgsea(m2014.de,c2.cp,protein_coding,"m2014")
b2007.fgsea    <- run_fgsea(b2007.de,c2.cp,protein_coding,"b2007")
c2004.fgsea    <- run_fgsea(c2004.de,c2.cp,protein_coding,"c2004")
s2012.fgsea    <- run_fgsea(s2012.de,c2.cp,protein_coding,"s2012")
ra2009a.fgsea  <- run_fgsea(ra2009a.de,c2.cp,protein_coding,"ra2009a")
ra2009b.fgsea  <- run_fgsea(ra2009b.de,c2.cp,protein_coding,"ra2009b")
n2010.fgsea    <- run_fgsea(n2010.de,c2.cp,protein_coding,"n2010")



sel <- m2014.de[m2014.de$logFC > 0 & m2014.de$P.Value <= 0.01,]$Gene_name
up1 <- unlist(strsplit(as.character(sel),',',fixed = T))
sel <- m2014.de[m2014.de$logFC < 0 & m2014.de$P.Value <= 0.01,]$Gene_name
dn1 <- unlist(strsplit(as.character(sel),',',fixed = T))
cat(sprintf("Up-regulated: %d, down-regulated: %d\n",length(up1),length(dn1)))

sel <- b2007.de[b2007.de$logFC > 0 & b2007.de$P.Value <= 0.01,]$Gene_name
up2 <- unlist(strsplit(as.character(sel),',',fixed = T))
sel <- b2007.de[b2007.de$logFC < 0 & b2007.de$P.Value <= 0.01,]$Gene_name
dn2 <- unlist(strsplit(as.character(sel),',',fixed = T))
cat(sprintf("Up-regulated: %d, down-regulated: %d\n",length(up2),length(dn2)))

sel <- c2004.de[c2004.de$logFC > 0 & c2004.de$P.Value <= 0.01,]$Gene_name
up3 <- unlist(strsplit(as.character(sel),',',fixed = T))
sel <- c2004.de[c2004.de$logFC < 0 & c2004.de$P.Value <= 0.01,]$Gene_name
dn3 <- unlist(strsplit(as.character(sel),',',fixed = T))
cat(sprintf("Up-regulated: %d, down-regulated: %d\n",length(up3),length(dn3)))

sel <- s2012.de[s2012.de$logFC > 0 & s2012.de$P.Value <= 0.01,]$Gene_name
up4 <- unlist(strsplit(as.character(sel),',',fixed = T))
sel <- s2012.de[s2012.de$logFC < 0 & s2012.de$P.Value <= 0.01,]$Gene_name
dn4 <- unlist(strsplit(as.character(sel),',',fixed = T))
cat(sprintf("Up-regulated: %d, down-regulated: %d\n",length(up4),length(dn4)))

sel <- ra2009a.de[ra2009a.de$logFC > 0 & ra2009a.de$P.Value <= 0.01,]$Gene_name
up5 <- unlist(strsplit(as.character(sel),',',fixed = T))
sel <- ra2009a.de[ra2009a.de$logFC < 0 & ra2009a.de$P.Value <= 0.01,]$Gene_name
dn5 <- unlist(strsplit(as.character(sel),',',fixed = T))
cat(sprintf("Up-regulated: %d, down-regulated: %d\n",length(up5),length(dn5)))

sel <- ra2009b.de[ra2009b.de$logFC > 0 & ra2009b.de$P.Value <= 0.01,]$Gene_name
up6 <- unlist(strsplit(as.character(sel),',',fixed = T))
sel <- ra2009b.de[ra2009b.de$logFC < 0 & ra2009b.de$P.Value <= 0.01,]$Gene_name
dn6 <- unlist(strsplit(as.character(sel),',',fixed = T))
cat(sprintf("Up-regulated: %d, down-regulated: %d\n",length(up6),length(dn6)))

sel <- n2010.de[n2010.de$logFC > 0 & n2010.de$P.Value <= 0.01,]$Gene_name
up7 <- unlist(strsplit(as.character(sel),',',fixed = T))
sel <- n2010.de[n2010.de$logFC < 0 & n2010.de$P.Value <= 0.01,]$Gene_name
dn7 <- unlist(strsplit(as.character(sel),',',fixed = T))
cat(sprintf("Up-regulated: %d, down-regulated: %d\n",length(up7),length(dn7)))

all_up <- unique(c(up1,up2,up3,up4,up5,up6,up7))
all_dn <- unique(c(dn1,dn2,dn3,dn4,dn5,dn6,dn6))

unfilt_up <- res_fann[res_fann$log2FoldChange > 0 & res_fann$padj <= 0.1,]
filt_up <- unfilt_up[! unfilt_up$Symbol %in% all_up,]
unfilt_dn <- res_fann[res_fann$log2FoldChange < 0 & res_fann$padj <= 0.1,]
filt_dn <- unfilt_dn[! unfilt_dn$Symbol %in% all_dn,]

write.table(filt_up,"filt_up_v1.tsv",quote = F,row.names = F,sep = "\t")
write.table(filt_dn,"filt_dn_v1.tsv",quote = F,row.names = F,sep = "\t")
write.table(unfilt_up,"unfilt_up_v1.tsv",quote = F,row.names = F,sep = "\t")
write.table(unfilt_dn,"unfilt_dn_v1.tsv",quote = F,row.names = F,sep = "\t")

m2014.de[m2014.de$Gene_name == "TLR9",]
n2010.de[n2010$Gene_name == "TLR9",]
n2010.de[n2010.de$Gene_name == "TLR9",]
ra2009a[ra2009a$Gene_name == "TLR9",]
b2007.de[b2007.de$Gene_name == "TLR9",]
c2004.de[c2004.de$Gene_name == "TLR9",]
s2012.de[s2012.de$Gene_name == "TLR9",]
ra2009a.de[ra2009a.de$Gene_name == "TLR9",]
ra2009a.de[ra2009a.de$Gene_name == "FOS",]
ra2009a.de[ra2009a.de$Gene_name == "CASP1",]
ra2009a.de[ra2009a.de$Gene_name == "NLRP3",]
ra2009b.de[ra2009b.de$Gene_name == "NLRP3",]

m2014.de[m2014.de$Gene_name == "HIF1A",]
n2010.de[n2010.de$Gene_name == "HIF1A",]
ra2009a.de[ra2009a.de$Gene_name == "HIF1A",]
ra2009b.de[ra2009b.de$Gene_name == "HIF1A",]
b2007.de[b2007.de$Gene_name == "HIF1A",]
c2004.de[c2004.de$Gene_name == "HIF1A",]
s2012.de[s2012.de$Gene_name == "HIF1A",]

m2014.de[m2014.de$Gene_name == "PHOSPHO1",]
n2010.de[n2010.de$Gene_name == "PHOSPHO1",]
ra2009a.de[ra2009a.de$Gene_name == "PHOSPHO1",]
ra2009b.de[ra2009b.de$Gene_name == "PHOSPHO1",]
b2007.de[b2007.de$Gene_name == "PHOSPHO1",]
c2004.de[c2004.de$Gene_name == "PHOSPHO1",]
s2012.de[s2012.de$Gene_name == "PHOSPHO1",]
