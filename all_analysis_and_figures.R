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
library(patchwork)
library(ggplotify)
library(cowplot)


setwd("~/tmp")
source("~/tmp/functions.R") 

################################## PART1. RNA-seq differential expression analysis #########################################

## PAR genes were removed from RSEM reference 
rsem_exp      <- read.table("rsem.genes.counts.tsv", header = T, row.names = 1, check.names = F)
rsem_tpm      <- read.table("rsem.genes.TPM.tsv", header = T, row.names = 1, check.names = F)

dim(rsem_exp)
head(rsem_exp)
head(rsem_tpm)

cond          <- read.table("Conditions.txt", row.names = 1, header = T, check.names = F,stringsAsFactors = T)
ann           <- rsem_exp[,c(1:2)]
exp           <- rsem_exp[,c(3:16)]
colSums(exp[, c(1:14)])

## make 18k table for visualization, do PCA on the rlog-transformed matrix, before and after donor correction with ComBat

exp_sorted      <- exp[order(rowMeans(exp), decreasing = T),]
exp_top18k      <- head(exp_sorted, 18000)
dds_top18k      <- DESeqDataSetFromMatrix(countData=round(exp_top18k,0),colData=cond,design = ~ donor + exercise)
mod18k          <- model.matrix(~ exercise, cond)
rlog_top18k     <- assay(rlog(dds_top18k,blind = F,fitType = "mean")) 
cbat_top18k     <- ComBat(rlog_top18k, batch = cond$donor, mod = mod18k)
boxplot(rlog_top18k)
boxplot(cbat_top18k)

p1 <- pca12t(rlog_top18k,cond,"exercise","donor","PCA before donor effect correction")
p2 <- pca12t(cbat_top18k,cond,"exercise","donor","PCA after donor effect correction") 

## Let's save annotated, ComBat-corrected table for heatmap visualizations

top18k_ann      <- merge(ann, cbat_top18k, by = 0)
colnames(top18k_ann)[1] <- c("Ensembl")
write.table(top18k_ann,"Skaters_top18k_cbat.tsv",quote = F,sep = "\t",row.names = F)

############# now let's do DESeq2 analysis and make some plots for diagnostics and illustration

dds_full      <- DESeqDataSetFromMatrix(countData=round(exp,0),colData = cond,design = ~ donor + exercise)
dds_full$exercise <- relevel(dds_full$exercise, ref = "before") ## all our log2FC are = (after/before)
dds_full      <- DESeq(dds_full)
res_full      <- results(dds_full, contrast=c("exercise","after","before"))
res_df        <- as.data.frame(res_full)
res_ann       <- merge(ann,res_df,by = 0)
res_ann$lfcSE <- NULL
colnames(res_ann)[1] <- c("Ensembl")
write.table(res_ann,"Skaters_full_diffexp.tsv",quote = F,sep = "\t",row.names = F)

## make MA plot
res_fann <- res_ann[complete.cases(res_ann),]
## DESeq2::plotMA(res_full, main = "MA plot, before vs. after exercise",ylim = c(-4, 4),
## colLine = "black",colSig = "red",colNonSig = "gray80")

res_ma <- as.data.frame(res_full[res_full$baseMean > 0,])
res_ma$regulation <- "not sig"
res_ma$regulation[res_ma$padj <= 0.1 & res_ma$log2FoldChange > 0] <- "up"
res_ma$regulation[res_ma$padj <= 0.1 & res_ma$log2FoldChange < 0] <- "dn"
table(res_ma$regulation)

p3 <- ggplot(res_ma,aes(x = log2(baseMean),y = log2FoldChange,color = regulation)) + geom_point(size = 0.4) + 
  scale_color_manual(values = c("not sig" = "#999999", "up" = "#E86850", "dn" = "#587498")) + 
  theme_bw() + ggtitle("MA plot, FDR 10% sign cutoff") + theme(panel.grid.major = element_line(size = 0.25)) +
  ylab("log fold change")


## annotate diffexp results with TPM, make filtered sets

rsem_tpm$Mean <- rowMeans(rsem_tpm[,3:16])
mean_tpm <- rsem_tpm[,17,drop = F]
colnames(mean_tpm) <- "Mean_TPM"
res_tpm <- merge(res_fann,mean_tpm,by.x = "Ensembl",by.y = "row.names")
res_tpm$significant <- "Not sig"
res_tpm$significant[res_tpm$padj <= 0.1] <- "FDR < 0.1"
## make Volcano plot with dot size proportional to expression
res_label <- subset(res_tpm, abs(log2FoldChange) > 2 | padj < 1e-8)
res_label <- res_label[! grepl("-",res_label$Symbol),]

res_tpm$regulation <- "not sig"
res_tpm$regulation[res_tpm$log2FoldChange > 0 & res_tpm$padj <= 0.1] <- "up"
res_tpm$regulation[res_tpm$log2FoldChange < 0 & res_tpm$padj <= 0.1] <- "dn"

p4 <- ggplot(res_tpm, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulation, size = log10(Mean_TPM+1)), stroke = 0) + 
  geom_text_repel(data = res_label, aes(label = Symbol), size = 3) + 
  scale_size_continuous(range = c(0.2,3)) + 
  scale_x_continuous(limits = c(-3.5, 3.5)) + 
  scale_color_manual(values = c("not sig" = "#999999", "up" = "#E86850", "dn" = "#587498")) + 
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

de_genes$cutoff <- c("0+", "10+", "100+", "1000+")
de_melt <- reshape2::melt(de_genes)
colnames(de_melt)[2] <- "regulation" 

p5 <- ggplot(de_melt, aes(x = cutoff, y = value, fill = regulation)) +
  geom_bar(stat="identity", position=position_dodge(0.4), width = 0.3)  + 
  scale_fill_manual(values = c("up" = "#E86850", "dn" = "#587498")) + 
  xlab("Mean expression cutoff (TPM)") + 
  ylab("Differentially expressed genes") + 
  ggtitle("DE genes, before vs. after exercise") + theme_bw() + 
  theme(panel.grid.major = element_line(size = 0.25)) 

## patchwork figure 1
p1235 <- (p1 + p2)/(p3 + p5)
p1235 | p4

## Save filtered differentially expressed gene sets

write.table(res_tpm[res_tpm$padj <= 0.1 & res_tpm$Gene_type == "protein_coding",],
            "Skaters_all_de.tsv",quote = F,sep = "\t",row.names = F)
write.table(res_tpm[res_tpm$padj <= 0.1 & res_tpm$Mean_TPM >= 10 & res_tpm$Gene_type == "protein_coding",],
            "Skaters_all_de_tpm10.tsv",quote = F,sep = "\t",row.names = F)
write.table(res_tpm[res_tpm$padj <= 0.1 & res_tpm$Mean_TPM >= 100 & res_tpm$Gene_type == "protein_coding",],
            "Skaters_all_de_tpm100.tsv",quote = F,sep = "\t",row.names = F)

################################## PART2. Pathway enrichment with enrichR and fGSEA #######################################

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

enr_up_h <- as.data.frame(enr_up_h)
enr_up_h <- enr_up_h[1:10,]
enr_up_h$fRatio <- sapply(enr_up_h$GeneRatio, function(x) eval(parse(text = x)))
colnames(enr_up_h)[9] <- "Overlap"
enr_up_h$ID <- gsub("_"," ",enr_up_h$ID)

enr_up_cp <- as.data.frame(enr_up_cp)
enr_up_cp <- enr_up_cp[1:10,]
enr_up_cp$fRatio <- sapply(enr_up_cp$GeneRatio, function(x) eval(parse(text = x)))
colnames(enr_up_cp)[9] <- "Overlap"
enr_up_cp$ID <- gsub("_"," ",enr_up_cp$ID)
enr_up_cp$ID[4] <- "REACTOME PLATELET ACTIVATION ..." ## too long of a name messes up the plot

enr_dn_cp <- as.data.frame(enr_dn_cp)
enr_dn_cp <- enr_dn_cp[1:10,]
enr_dn_cp$fRatio <- sapply(enr_dn_cp$GeneRatio, function(x) eval(parse(text = x)))
colnames(enr_dn_cp)[9] <- "Overlap"
enr_dn_cp$ID <- gsub("_"," ",enr_dn_cp$ID)
enr_dn_cp$ID[2] <- "REACTOME RRNA MODIFICATION ..."

q1 <- ggplot(enr_up_h,aes(x = fRatio,y = factor(ID,levels = ID[order(fRatio)]),size = Overlap,color = -log(qvalue))) + 
  geom_point() + scale_color_gradient(low = "blue",high = "red") + 
  xlab("Gene ratio") + 
  ggtitle("Gene overlap, up-regulated genes vs. MsigDB H") + theme_bw() + theme(plot.title = element_text(hjust = 1))
q2 <- ggplot(enr_up_cp,aes(x = fRatio,y = factor(ID,levels = ID[order(fRatio)]),size = Overlap,color = -log(qvalue))) + 
  geom_point() + scale_color_gradient(low = "blue",high = "red") + 
  xlab("Gene ratio") + 
  ggtitle("Gene overlap, up-regulated genes vs. MsigDB C2:CP") + theme_bw() + theme(plot.title = element_text(hjust = 1))
q3 <- ggplot(enr_dn_cp,aes(x = fRatio,y = factor(ID,levels = ID[order(fRatio)]),size = Overlap,color = -log(qvalue))) + 
  geom_point() + scale_color_gradient(low = "blue",high = "red") + 
  xlab("Gene ratio") + 
  ggtitle("Gene overlap, down-regulated genes vs. MsigDB C2:CP") + theme_bw() + theme(plot.title = element_text(hjust = 1))

#q1 <- dotplot(enr_up_h, showCategory = 10, title = "Gene overlap, up-regulated genes vs. MsigDB H")
#q2 <- dotplot(enr_up_cp, showCategory = 10,title = "Gene overlap, up-regulated genes vs. MsigDB C2:CP")
#q3 <- dotplot(enr_dn_cp, showCategory = 10,title = "Gene overlap, down-regulated genes vs. MsigDB C2:CP")

## run similar enrichments using fGSEA
rnk            <- aggregate(stat ~ Symbol,data = res_tpm[res_tpm$Gene_type == "protein_coding",], 
                            function(x) ifelse(mean(x)>0,max(x),min(x)))
rnk            <- setNames(rnk$stat,toupper(rnk$Symbol))
length(rnk)

gsea.h         <- fgsea(h.all,rnk,minSize = 10,maxSize = 500,eps = 0.0,nproc = 4)
gsea.cp        <- fgsea(c2.cp,rnk,minSize = 10,maxSize = 500,eps = 0.0,nproc = 4)

gsea_up_h <- gsea.h[gsea.h$NES > 0 & gsea.h$padj < 0.05,]
gsea_up_h <- gsea_up_h[order(gsea_up_h$padj),]
gsea_up_h <- gsea_up_h[1:10,]
gsea_up_h$Overlap <- lengths(gsea_up_h$leadingEdge)
gsea_up_h$pathway <- gsub("_"," ",gsea_up_h$pathway)

gsea_dn_h <- gsea.h[gsea.h$NES < 0 & gsea.h$padj < 0.05,] 
gsea_dn_h <- gsea_up_h[order(gsea_up_h$padj),] ## very similar situation, only MYC v2 is very significant
gsea_up_cp <- gsea.cp[gsea.cp$NES > 0 & gsea.cp$padj < 0.05,]
gsea_up_cp <- gsea_up_cp[order(gsea_up_cp$padj),]
gsea_up_cp <- gsea_up_cp[1:10,]
gsea_up_cp$Overlap <- lengths(gsea_up_cp$leadingEdge)
gsea_up_cp$pathway <- gsub("_"," ",gsea_up_cp$pathway)
gsea_up_cp$pathway[8] <- "REACTOME DISEASES OF SIGNAL ..."
gsea_up_cp$pathway[3] <- "REACTOME PLATELET ACTIVATION ..."


gsea_dn_cp <- gsea.cp[gsea.cp$NES < 0 & gsea.cp$padj < 0.05,] 
gsea_dn_cp <- gsea_dn_cp[order(gsea_dn_cp$padj),] 
gsea_dn_cp <- gsea_dn_cp[1:10,]
gsea_dn_cp$Overlap <- lengths(gsea_dn_cp$leadingEdge)
gsea_dn_cp$pathway <- gsub("_"," ",gsea_dn_cp$pathway)
gsea_dn_cp$pathway[9] <- "REACTOME SRP DEPENDENT ..."
gsea_dn_cp$pathway[8] <- "REACTOME RRNA MODIFICATION ..."


## let's plot same enrichments as above but for GSEA
q4 <- ggplot(gsea_up_h,aes(x = NES,y = factor(pathway, levels = gsea_up_h$pathway[order(gsea_up_h$NES)]),
                           size = Overlap,color = -log(padj))) + 
  geom_point() + scale_color_gradient(low = "blue",high = "red") + 
  xlab("Normalized enrichment score (NES)") + 
  ggtitle("GSEA enrichment, up-regulated pathways from MsigDB H") + theme_bw() + theme(plot.title = element_text(hjust = 1))
q5 <- ggplot(gsea_up_cp,aes(x = NES,y = factor(pathway, levels = gsea_up_cp$pathway[order(gsea_up_cp$NES)]),
                            size = Overlap,color = -log(padj))) + 
  geom_point() + scale_color_gradient(low = "blue",high = "red") + 
  xlab("Normalized enrichment score (NES)") + 
  ggtitle("GSEA enrichment, up-regulated pathways from MsigDB CP") + theme_bw() + theme(plot.title = element_text(hjust = 1))
q6 <- ggplot(gsea_dn_cp,aes(x = -NES,y = factor(pathway, levels = gsea_dn_cp$pathway[order(-gsea_dn_cp$NES)]),
                            size = Overlap,color = -log(padj))) + 
  geom_point() + scale_color_gradient(low = "blue",high = "red") + 
  xlab("Negative normalized enrichment score (-NES)") + 
  ggtitle("GSEA enrichment, down-regulated pathways from MsigDB CP") + theme_bw() + theme(plot.title = element_text(hjust = 1))

## patchwork figure 2; export 8 x 20 pdf
(q1 + q2 + q3)/(q4 + q5 + q6) & theme(axis.title.y = element_blank(),plot.margin = margin(30, 30, 30, 30))

plot_grid(q1, q2, q3, q4, q5, q6, nrow=2)
################################## PART3. Puting our dataset in the context of other (microarray) datasets ########################


## use GEOquery package to get the selected datasets
m2014.gse     <- getGEO("GSE51216",GSEMatrix = T)[[1]] ## GPL6480, log2-tr, normalized
t2013.gse     <- getGEO("GSE46075",GSEMatrix = T)[[1]] ## GPL6244, not log2, normalized
b2007.gse     <- getGEO("GSE3606",GSEMatrix = T)[[1]] ## GPL571, not log2, normalized
c2004.gse     <- getGEO("GSE1140",GSEMatrix = T)[[1]] ## GPL96, not log2, normalized
s2012.gse     <- getGEO("GSE28498",GSEMatrix = T)[[1]] ## GPL6244, not log2, normalized (t2 & t3 have clear batch)
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
          "GSM706111","GSM706112","GSM706113","GSM706114","GSM706115","GSM706116","GSM706117","GSM706118","GSM706119",
          "GSM706120","GSM706121")
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

## custom annotations are used with updated gene names etc (all are in /data folder)
gpl96.ann    <- read.table("GPL96.ann.tsv",header = T, row.names = 1)
gpl570.ann   <- read.table("GPL570.ann.tsv",header = T, row.names = 1)
gpl571.ann   <- read.table("GPL571.ann.tsv",header = T, row.names = 1)
gpl6244.ann  <- read.table("GPL6244.ann.tsv",header = T, row.names = 1)
gpl6480.ann  <- read.table("GPL6480.ann.tsv",header = T, row.names = 1)
gpl4133.ann  <- read.table("GPL4133.ann.tsv",header = T, row.names = 1)

## limma_donor_correct is available in functions.R
##  most of these look amazingly consistent. 
m2014.de     <- limma_donor_correct(m2014.gse,m2014.cond,gpl6480.ann,log2transformed = T,normalized = T,"Mukherjee 2014, GSE51216")
t2013.de     <- limma_donor_correct(t2013.gse,t2013.cond,gpl6244.ann,log2transformed = F,normalized = F,"Tonevitsky 2013, GSE46075") 
b2007.de     <- limma_donor_correct(b2007.gse,b2007.cond,gpl571.ann,log2transformed = F,normalized = T,"Buttner 2007, GSE3606") 
c2004.de     <- limma_donor_correct(c2004.gse,c2004.cond,gpl96.ann,log2transformed = F,normalized = T,"Connoly 2004, GSE1140")
s2012.de     <- limma_donor_correct(s2012.gse,s2012.cond,gpl6244.ann,log2transformed = F,normalized = T,"Sakharov 2012, GSE28498")
ra2009a.de   <- limma_donor_correct(ra2009a.gse,ra2009a.cond,gpl570.ann,log2transformed = F,normalized = T,"Radom-Azik 2009a, GSE14642")
ra2009b.de   <- limma_donor_correct(ra2009b.gse,ra2009b.cond,gpl570.ann,log2transformed = F,normalized = T,"Radom-Azik 2009b, GSE11761")
n2010.de     <- limma_donor_correct(n2010.gse,n2010.cond,gpl4133.ann,log2transformed = F,normalized = F,"Nakamura 2010, GSE18966")


## now, we take all up-regulated genes (relaxed criteria, non-adjusted pval <= 0.01) and subtract them from the RNA-seq DE genes
## we don't use Tonevitsky 2013 since it looks strange (as if samples were mixed up)
de_all <- do.call("rbind", list(m2014.de,b2007.de,c2004.de,s2012.de,ra2009a.de,ra2009b.de,n2010.de))
dim(de_all)

array_up <- de_all[de_all$logFC > 0 & de_all$P.Value <= 0.01,]$Gene_name
array_up <- unique(unlist(strsplit(as.character(array_up),',',fixed = T)))
length(array_up)
array_dn <- de_all[de_all$logFC < 0 & de_all$P.Value <= 0.01,]$Gene_name
array_dn <- unique(unlist(strsplit(as.character(array_dn),',',fixed = T)))
length(array_dn)

## remove all genes that were up or down-regulated in the 7 microarray studies
rnaseq_up <- res_tpm[res_tpm$log2FoldChange > 0 & res_tpm$padj <= 0.1,]
filt_up <- rnaseq_up[! rnaseq_up$Symbol %in% array_up,]
rnaseq_dn <- res_tpm[res_tpm$log2FoldChange < 0 & res_tpm$padj <= 0.1,]
filt_dn <- rnaseq_dn[! rnaseq_dn$Symbol %in% array_dn,]

## These are all the genes that are differentially expressed in RNA-seq 
## and not detected changing in the same direction in any of the evaluated microarray datasets.
filtered <- rbind(filt_up,filt_dn)
write.table(filtered,"Array_filtered_de.tsv",quote = F,row.names = F,sep = "\t")

################################## PART4. Single-cell RNA-seq datasets used to generate markers ########################

library(Seurat)
library(ggplot2)
library(dplyr)
library(DropletUtils)
library(hdf5r)

## Very few datasets have erythrocytes and/or whole blood. We'll use this one:
## this is dataset https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149938

download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149938/suppl/GSE149938_umi_matrix.csv.gz",
              destfile = "GSE149938_umi_matrix.csv.gz")
umi_gz <- gzfile("GSE149938_umi_matrix.csv.gz",'rt')  
umi <- read.csv(umi_gz,check.names = F,quote = "")

wb <- CreateSeuratObject(t(umi),project = "whole_blood") ## works nicely! check out orig.ident
wb[["percent.mt"]] <- PercentageFeatureSet(wb, pattern = "^MT-") ## apparently no mito genes? 
wb[["percent.rbp"]] <- PercentageFeatureSet(wb, pattern = "^RP[SL]")
VlnPlot(wb, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

## process all 
wb <- NormalizeData(wb)
wb <- FindVariableFeatures(wb, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(wb)
wb <- ScaleData(wb, features = all.genes)
wb <- RunPCA(wb, features = VariableFeatures(object = wb))
DimPlot(wb, reduction = "pca")
ElbowPlot(wb)

## check top genes. As expected, it's mostly globins and ribosomal proteins
feat <- wb[['RNA']]@meta.features
feat <- feat[order(feat$vst.mean,decreasing = T),]
head(feat,n = 50)

## clustering at low resolution is great to define bigger cell populations
wb <- RunUMAP(wb, dims = 1:20)
wb <- FindNeighbors(wb, dims = 1:20)
wb <- FindClusters(wb,resolution = 0.1) 
DimPlot(wb,label.size = 6,repel = T,label = T)

## Let's see the original annotations. 
wb <- SetIdent(wb,value = "orig.ident")
DimPlot(wb,label.size = 6,repel = T,label = T)
wb <- SetIdent(wb,value = "seurat_clusters")


## putative neutrophil markers
FeaturePlot(wb,c("CEACAM8","CAMP","LCN2","LTF"))
## putative erythrocyte markers
FeaturePlot(wb,c("HBA1","HBA2","HBB","HBD"))

## We have confrimed that cluster 0 is erythrocytes, cluster 6 is neutrophils. 
## Let's define the markers for these cell types. 

wb.markers <- FindAllMarkers(wb, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

#ery <- wb.markers[wb.markers$cluster == "0",]
#ery <- ery[! grepl("^HB",ery$gene),] ## lose the globins since we depleted them
#ery$cell <- "erythrocyte"
#ery <- ery[1:20,]

#neutro <- wb.markers[wb.markers$cluster == "6",]
#neutro$cell <- "neutrophil"
#neutro <- neutro[1:20,]

### Let's process 10X PBMC10k dataset the same way and get all other markers 
download.file("https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_filtered_feature_bc_matrix.h5",
              destfile = "pbmc10k_filt.h5")
filt.matrix <- Read10X_h5("pbmc10k_filt.h5",use.names = T)
pbmc        <- CreateSeuratObject(counts = filt.matrix)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 15) ## crude filtering for doublets
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc    <- NormalizeData(pbmc)
pbmc    <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

## Again, we want courser clustering - overall B/T/Monocyte/NK/platelet 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:15)
pbmc <- FindNeighbors(pbmc, dims = 1:15, k.param = 15)
pbmc <- FindClusters(pbmc, resolution = 0.2)

DimPlot(pbmc,label.size = 6,repel = T,label = T)
FeaturePlot(pbmc,c("LILRA4"))
FeaturePlot(pbmc,c("FCER1A"))
FeaturePlot(pbmc,c("PPBP")) ## got to split clusters, these are DCs and platelets together
FeaturePlot(pbmc,c("NKG7"))

plat_barcodes <- colnames(subset(x = pbmc, subset = PPBP > 2))
#miel_barcodes <- colnames(subset(x = pbmc, subset = FCER1A > 2))
pbmc@meta.data$seurat_clusters <- as.character(pbmc@meta.data$seurat_clusters)
pbmc@meta.data$seurat_clusters[pbmc@meta.data$seurat_cluster == '7'] <- '5'
pbmc@meta.data$seurat_clusters[rownames(pbmc@meta.data) %in% plat_barcodes] <- '7'
#pbmc@meta.data$seurat_clusters[rownames(pbmc@meta.data) %in% miel_barcodes] <- '9'
table(pbmc@meta.data$seurat_clusters)
pbmc@meta.data$seurat_clusters <- as.factor(pbmc@meta.data$seurat_clusters)
pbmc <- SetIdent(pbmc, value = "seurat_clusters")
DimPlot(pbmc,label.size = 6,repel = T,label = T)

## Now we're ready to define the markers for 8 populations. Cluster identification can be done using 
## markers listed in Seurat vignette: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
pbmc.markers <- FindAllMarkers(pbmc, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers <- pbmc.markers[! grepl("^RP[SL]",pbmc.markers$gene),]
pbmc.markers <- pbmc.markers[! grepl("\\.",pbmc.markers$gene),]
pbmc.markers <- pbmc.markers[! grepl("^LINC",pbmc.markers$gene),]
# Examine the number of clusters per marker - expecting one cluster per gene
table(table(pbmc.markers$gene))

# Exploring known markers
FeaturePlot(pbmc,c("CD14"))
FeaturePlot(pbmc,c("IL7R"))
FeaturePlot(pbmc,c("GATA3"))
FeaturePlot(pbmc,c("CD8A"))
FeaturePlot(pbmc,c("IGLC2"))
FeaturePlot(pbmc,c("TCF7L2"))
FeaturePlot(pbmc,c("GZMH"))
DimPlot(pbmc,label.size = 6,repel = T,label = T)

# Combine WB and pbmc
wb.markers$dataset = 'wb'
pbmc.markers$dataset = 'pbmc'
all.markers = rbind(wb.markers[wb.markers$cluster %in% c('0', '6'), ], 
                    pbmc.markers)
table(all.markers$cluster)
all.markers$cluster = ifelse(all.markers$dataset == 'wb',
                  ifelse(all.markers$cluster == '0', '10', '11'),
                  as.character(all.markers$cluster))
table(all.markers$cluster)
head(all.markers)
table(table(all.markers$gene))
all.uq.markers = all.markers[all.markers$gene %in% 
          names(table(all.markers$gene))[table(all.markers$gene) == 1], ]
table(all.uq.markers$cluster)

ery <- all.uq.markers[all.uq.markers$cluster == "10",]
ery <- ery[! grepl("^HB",ery$gene),] ## lose the globins since we depleted them
ery$cell <- "erythrocyte"
ery <- ery[1:20,]

neutro <- all.uq.markers[all.uq.markers$cluster == "11",]
neutro$cell <- "neutrophil"
neutro <- neutro[1:20,]

cd14_mono <- all.uq.markers[all.uq.markers$cluster == "0",]
cd14_mono$cell <- "CD14_monocyte"
cd14_mono <- cd14_mono[1:20,]

cd4_naive <- all.uq.markers[all.uq.markers$cluster == "1",]
cd4_naive$cell <- "CD4_naive_Tcell"
cd4_naive <- cd4_naive[1:20,]

cd4_memory <- all.uq.markers[all.uq.markers$cluster == "2",]
cd4_memory$cell <- "CD4_memory_Tcell"
cd4_memory <- cd4_memory[1:20,]

cd8_tcell <- all.uq.markers[all.uq.markers$cluster == "3",]
cd8_tcell$cell <- "CD8_Tcell"
cd8_tcell <- cd8_tcell[1:20,]

nk <- all.uq.markers[all.uq.markers$cluster == "5",]
nk$cell <- "natural_killer"
nk <- nk[1:20,]

bcell <- all.uq.markers[all.uq.markers$cluster == "4",]
bcell$cell <- "Bcell"
bcell <- bcell[1:20,]

cd16_mono <- all.uq.markers[all.uq.markers$cluster == "6",]
cd16_mono$cell <- "CD16_monocyte"
cd16_mono <- cd16_mono[1:20,]

pldc <- all.uq.markers[all.uq.markers$cluster == "8",]
pldc$cell <- "plasmacytoid_dendritic"
pldc <- pldc[1:20,]

mydc <- all.uq.markers[all.uq.markers$cluster == "9",]
mydc$cell <- "myeloid_dendritic"
mydc <- mydc[1:20,]

plat <- all.uq.markers[all.uq.markers$cluster == "7",]
plat$cell <- "platelet"
plat <- plat[1:20,]

all_top20_markers <- rbind(ery,neutro,cd14_mono,cd16_mono,cd4_naive,cd4_memory,cd8_tcell,nk,bcell,mydc,pldc,plat)
## replace some of the outdated gene names 
all_top20_markers$gene[! all_top20_markers$gene %in% res_ann$Symbol]
## "FAM101B" -> "RFLNB"; "CAVIN2" -> "SDPR"
all_top20_markers$gene[all_top20_markers$gene == "FAM101B"] <- "RFLNB"
all_top20_markers$gene[all_top20_markers$gene == "CAVIN2"]  <- "SDPR"

## make the table of markers used
write.table(all_top20_markers,"Top20_markers_scRNAseq.tsv",quote = F,row.names = F,sep = "\t")

### now, let's make a 


################################## PART5. Identify changes of cell populations in bulk RNA-seq ########################

library(pheatmap)
library(ComplexHeatmap)

## Let's take expression matrix after donor correction, and find all the defined scRNAseq cell-type specific markers. 
cbat_select <- top18k_ann[top18k_ann$Symbol %in% all_top20_markers$gene,]
dim(cbat_select)
cbat_markers <- merge(all_top20_markers[,c(7,9)],cbat_select,by.x = "gene",by.y = "Symbol")
cbat_markers$Ensembl <- NULL
cbat_markers$Gene_type <- NULL
markers_melt <- melt(cbat_markers)
markers_melt$exercise <- gsub("\\d","",markers_melt$variable,perl = T)
markers_melt$exercise <- factor(markers_melt$exercise,levels = c("before","after"))
markers_melt$variable <- NULL

markers_melt$cell = factor(markers_melt$cell, levels=c('erythrocyte',
        'CD14_monocyte', 'neutrophil', 'platelet', 'CD4_memory_Tcell',
        'CD8_Tcell', 'Bcell', 'natural_killer', 'plasmacytoid_dendritic',
        'CD4_naive_Tcell', 'CD16_monocyte', 'myeloid_dendritic'))

stat_data <- data.frame(cell=levels(markers_melt$cell),
                        wilcox_down=NA, wilcox_up=NA,
                        cell_index=1:12)
rownames(stat_data) = stat_data$cell

# Statistical comparison using paired Wilcoxon criterion
for (cell in unique(markers_melt$cell)) {
  expr_cell <- aggregate(value~gene+exercise, 
                         markers_melt[as.character(markers_melt$cell) ==
                                        as.character(cell), ], median)
  stat_data[cell, 'wilcox_down'] = wilcox.test(value~exercise, expr_cell, paired=T, alternative="greater")$p.value
  stat_data[cell, 'wilcox_up'] = wilcox.test(value~exercise, expr_cell, paired=T, alternative="less")$p.value
  print(paste0(cell, 
               " - p (down) = ", 
               wilcox.test(value~exercise, expr_cell, paired=T, alternative="greater")$p.value,
               ", p (up) = ", 
               wilcox.test(value~exercise, expr_cell, paired=T, alternative="less")$p.value))
}

stat_data$annotation = ifelse(stat_data$wilcox_up < 0.001 | stat_data$wilcox_down < 0.001, 
                              "***", ifelse(stat_data$wilcox_up < 0.01 | stat_data$wilcox_down < 0.01,
                                            "**", "-"))

library(ggsignif)
stat_data$x_start = stat_data$cell_index - 0.25
stat_data$x_end = stat_data$cell_index + 0.25
stat_data$y = 18

## boxplot of all markers, before and after exercise
r1_box <- ggplot(markers_melt,aes(x=cell,y=value)) + 
  geom_violin(aes(fill=exercise), scale='width') + 
  geom_signif(stat="identity", 
              data=stat_data, 
              aes(x=x_start, xend=x_end, 
                  y=y, yend=y, annotation=annotation, group=cell_index)) +
  geom_boxplot(aes(fill=exercise), width=0.2, outlier.shape=NA, position=position_dodge(width=0.9)) +
  theme_bw() + scale_fill_manual(values = c("#587058","#FFD800")) + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  ylab("Top 20 single cell markers") + xlab("cell type") +
  scale_y_continuous(limits=c(0, 20))
r1_box

markers_melt$group <- paste0(markers_melt$cell,markers_melt$gene)
markers_melt$group <- as.factor(markers_melt$group)

##   geom_boxplot(width = 0.6,position=position_dodge(1)) + 
r1 <- ggplot(markers_melt,aes(x = cell,y = value,fill = exercise)) + theme_bw() + 
  geom_path(position=position_dodge(1),group = markers_melt$group, size = 0.1) + 
  geom_point(aes(color = exercise), colour="black",pch = 21,stroke = 0.2, position=position_dodge(1), size = 2) + 
  scale_color_manual(values = c("#587058","#FFD800")) + scale_fill_manual(values = c("#587058","#FFD800")) + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  ylab("Top 20 single cell markers") + xlab("cell type") + ggtitle("Expression changes of top 20 cell type specific markers")





## Let's now make a heatmap for all markers
cell_types <- unique(cbat_markers$cell)
hmap <- list()
sorted_cond <- cond[order(cond$exercise,decreasing = T),1,drop = F]

for (i in cell_types) {
  cat(sprintf("Processing cell type: %s\n",i))
  marker_exp <- cbat_markers[cbat_markers$cell == i,]
  rownames(marker_exp) <- marker_exp$gene
  marker_exp$gene <- NULL
  marker_exp <- marker_exp[,rownames(sorted_cond)]
  hmap[[i]] <- as.ggplot(pheatmap(marker_exp,cluster_rows = F,annotation_col = sorted_cond,cluster_cols = F,scale = "row",
                     labels_col = "", main = i,legend = F,annotation_legend = F,border_color = "white", 
                     annotation_colors = list(exercise = c(after = "#FFD800",before = "#587058"))))
}

r2 <- wrap_elements(grid::textGrob('Increased')) + hmap[["erythrocyte"]] + 
  hmap[["CD14_monocyte"]] + hmap[["neutrophil"]] + hmap[['platelet']] + plot_layout(ncol = 5, widths = c(1,3,3,3,3))
r3 <- wrap_elements(grid::textGrob('Decreased')) + hmap[["CD4_memory_Tcell"]] + 
  hmap[["CD8_Tcell"]] + hmap[["Bcell"]] + hmap[["natural_killer"]] + plot_layout(ncol = 5, widths = c(1,3,3,3,3))
r4 <- wrap_elements(grid::textGrob('Unchanged')) + hmap[["CD4_naive_Tcell"]] + 
  hmap[["CD16_monocyte"]] + hmap[["myeloid_dendritic"]] + hmap[["plasmacytoid_dendritic"]] + plot_layout(ncol = 5, widths = c(1,3,3,3,3))


r1 <- ggplot(markers_melt,aes(x = cell,y = value,fill = exercise)) + theme_bw() + 
  geom_path(position=position_dodge(1),group = markers_melt$group, size = 0.1) + 
  geom_point(aes(color = exercise), colour="black",pch = 21,stroke = 0.2, position=position_dodge(1), size = 2) + 
  scale_color_manual(values = c("#587058","#FFD800")) + scale_fill_manual(values = c("#587058","#FFD800")) + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  ylab("Top 20 single cell markers") + xlab("cell type") + ggtitle("Expression changes of top 20 cell type specific markers")


##   geom_boxplot(width = 0.6,position=position_dodge(1)) + 
cairo_pdf('Figure_3_revamp.pdf', width=18, height=15)
plot_grid(r1_box, r2, r3, r4, nrow=4, rel_heights = c(1, 0.8, 0.8, 0.8))
dev.off()

## Patchwork figure 3, 15 x 18. TODO: sort cells by unchanged/dec/inc. 
r1/r2/r3/r4


################################## PART6. Annotate and plot more "array-filtered" genes ########################

## This is a logical continuation of PART3, but with blackjack and single cell markers. 
dim(filtered)
head(filtered)


table(all_top20_markers$cell)
filt_ann <- filtered
filt_ann$cell <- "Unknown"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "1",]$gene] <- "CD4_naive_Tcell"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "2",]$gene] <- "CD4_memory_Tcell"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "3",]$gene] <- "CD8_Tcell"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "4",]$gene] <- "Bcell"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "5",]$gene] <- "natural_killer"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "6",]$gene] <- "CD16_monocyte"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "8",]$gene] <- "plasmacytoid_dendritic"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "9",]$gene] <- "myeloid_dendritic"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "7",]$gene] <- "platelet"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "0",]$gene] <- "CD14_monocyte"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "10",]$gene] <- "erythrocyte"
filt_ann$cell[filt_ann$Symbol %in% all.uq.markers[all.uq.markers$cluster == "11",]$gene] <- "neutrophil"


filt_ann$cell[filt_ann$Symbol %in% wb.markers[wb.markers$cluster == "0",]$gene] <- "erythrocyte"
filt_ann$cell[filt_ann$Symbol %in% wb.markers[wb.markers$cluster == "6",]$gene] <- "neutrophil"
table(filt_ann$cell)
write.table(filt_ann,"Full_filtered_celltype_de.tsv",quote = F,row.names = F,sep = "\t")

filt_label <- subset(filtered, Mean_TPM > 10 & Gene_type == "protein_coding" & abs(log2FoldChange) >= 0.5)


### We can make a volcano plot of the filtered genes, with all other DE genes in the background. 

res_array <- res_tpm[res_tpm$significant == "FDR < 0.1",]
res_array$filter <- "array-filtered"
res_array$filter[res_array$Symbol %in% filtered$Symbol[filtered$regulation == "up"]] <- "up (RNAseq)"
res_array$filter[res_array$Symbol %in% filtered$Symbol[filtered$regulation == "dn"]] <- "dn (RNAseq)"
table(res_array$filter)

filt_label2 <- subset(filtered, Mean_TPM > 5)
filt_label2 <- subset(filt_label2, padj < 1e-4 | abs(log2FoldChange) >= 0.5)
filt_label2 <- filt_label2[! grepl("-",filt_label2$Symbol),] ## don't want -AS* etc

write.table(filt_label2,"Array_filtered_celltype_de.tsv",quote = F,row.names = F,sep = "\t")


## select protein coding genes for heatmap vis
filt_label3 <- subset(filt_ann, Mean_TPM > 10 & abs(log2FoldChange) >= 0.5 & Gene_type == "protein_coding")


s1 <- ggplot(res_array, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(res_array, filter = "array-filtered"), aes(color = filter, size = log10(Mean_TPM+1)), stroke = 0) + 
  geom_point(data = subset(res_array, filter != "array-filtered"), aes(color = filter, size = log10(Mean_TPM+1)), stroke = 0) + 
  geom_text_repel(data = filt_label2, aes(label = Symbol), size = 3) + 
  scale_size_continuous(range = c(0.2,3)) + 
  scale_x_continuous(limits = c(-2, 2)) + 
  scale_color_manual(values = c("array-filtered" = "#ADD8E6", "up (RNAseq)" = "#E86850", "dn (RNAseq)" = "#587498")) + 
  xlab("log fold change") + 
  ylab("-log10 (adjusted p-value)") + 
  ggtitle("Genes uniquely up-/down-regulated in our dataset") + theme_bw() + 
  theme(panel.grid.major = element_line(size = 0.25))


## Now let's make annotated heatmap. First, get the 72 markers
marker_exp <- top18k_ann[top18k_ann$Symbol %in% filt_label3$Symbol,]
rownames(marker_exp) <- marker_exp$Symbol
marker_exp[,1:3] <- NULL
marker_exp <- marker_exp[,rownames(sorted_cond)]

sorted_rows <- filt_label3[,c(2,12,5,9)]
rownames(sorted_rows) <- sorted_rows$Symbol
sorted_rows$Symbol <- NULL

## Order by cell type and then by log fold change
sorted_rows <- sorted_rows[order(sorted_rows$cell,sorted_rows$log2FoldChange,decreasing = T),]
marker_exp <- marker_exp[rownames(sorted_rows),]

## Adjust dataframes a bit more 
sorted_tpms <- sorted_rows[,3,drop = F]
sorted_tpms$Symbol <- rownames(sorted_tpms)
sorted_rows[,2:3] <- NULL
colnames(sorted_rows) <- "celltype"

sorted_tpms$Exp <- "TPM > 10"
sorted_tpms$Exp[sorted_tpms$Mean_TPM >= 50] <- "TPM > 50"
sorted_tpms$Exp[sorted_tpms$Mean_TPM >= 100] <- "TPM > 100"
sorted_tpms$Exp[sorted_tpms$Mean_TPM >= 500] <- "TPM > 500"

table(sorted_tpms$Exp)
sorted_tpms$Exp <- factor(sorted_tpms$Exp,levels = c("TPM > 10","TPM > 50","TPM > 100","TPM > 500"))
sorted_tpms$Symbol <- factor(sorted_tpms$Symbol,levels = rev(sorted_tpms$Symbol)) 

s2 <- ggplot(sorted_tpms,aes(x = Symbol,y = log10(Mean_TPM),fill = Exp)) + geom_bar(stat = "identity",width = 0.9) + 
  theme_bw() + coord_flip() + scale_y_continuous(position = "right") + scale_y_reverse() + 
  scale_fill_manual(values = c("TPM > 10" = "#fef0d9", "TPM > 50" = "#fdcc8a", "TPM > 100" = "#fc8d59","TPM > 500" = "#d7301f")) + 
  theme(legend.position = "bottom",panel.grid.major.y = element_blank()) + 
  ggtitle("Array-filtered genes, TPM > 10, |logFC| > 0.5")

s3 <- as.ggplot(pheatmap(marker_exp,cluster_rows = F,annotation_col = sorted_cond,annotation_row = sorted_rows, 
                   cluster_cols = F,scale = "row",labels_col = "", fontsize_row = 9,
                   legend = F,annotation_legend = T,border_color = "white", gaps_row = c(27,34,50,56,57,58,59,62,71,72),
                   annotation_names_row = F, annotation_names_col = F,
                   annotation_colors = list(exercise = c(after = "#FFD800",before = "#587058"))))

markers_matrix <- as.matrix(marker_exp)
scaled_matrix = (markers_matrix - rowMeans(markers_matrix))/apply(markers_matrix, 1, sd)

column_ha = ComplexHeatmap::HeatmapAnnotation(exercise = sorted_cond$exercise,
                      simple_anno_size = unit(0.4, "cm"),
                      col = list(exercise = c(after = "#FFD800",before = "#587058")))

all_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
bar_color = sapply(log10(sorted_tpms$Mean_TPM - min(sorted_tpms$Mean_TPM))/
                   log10(max(sorted_tpms$Mean_TPM) - min(sorted_tpms$Mean_TPM)),
            function(x) all_colors[ceiling(x * 100)][1])
factor_colors = c("TPM > 10" = "#fef0d9", "TPM > 50" = "#fdcc8a", "TPM > 100" = "#fc8d59","TPM > 500" = "#d7301f")
bar_factor_color = sapply(as.character(sorted_tpms$Exp),
                  function(x) factor_colors[x])
names(bar_factor_color) <- NULL
row_colors = list("cell" = rev(c(RColorBrewer::brewer.pal(n=12, name="Set3"), '#00d65c')))
names(row_colors['cell'][[1]]) = unique(filt_ann$cell)

row_ha = ComplexHeatmap::rowAnnotation(cell = sorted_rows$celltype,
            logtpm = ComplexHeatmap::anno_barplot(log10(sorted_tpms$Mean_TPM),
              width = unit(2, "cm"), gp=gpar(fill=bar_factor_color)),
            simple_anno_size = unit(0.5, "cm"), col=row_colors)
chmap <- ComplexHeatmap::Heatmap(scaled_matrix, name = "mat", 
        top_annotation = column_ha, right_annotation = row_ha,
        cluster_rows = F, cluster_columns = F, 
        row_split = sorted_rows$celltype, show_column_names = F, row_title = NULL, 
        col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
        row_names_gp = gpar(fontsize = 8), rect_gp = gpar(col = "white", lwd = 1))
chmap
s3_ch <- as.ggplot(chmap)
  
## patchwork figure4, 10 x 20 exp
cairo_pdf('Figure_4_revamp.pdf', width=21, height=9)
plot_grid(s1, s3_ch, nrow=1, rel_widths = c(1.1, 1))
dev.off()

s23 <- s2 + s3 + plot_layout(widths = c(1,5))
s1 + s23

# Computing the enrichment of markers in DEGs and unique DEGs
degs_rnaseq <- rbind(rnaseq_up, rnaseq_dn)
high_degs_rnaseq <- subset(degs_rnaseq, Mean_TPM > 10 & Gene_type == "protein_coding" & abs(log2FoldChange) >= 0.5)

marker_enrich <- data.frame(cell = unique(all_top20_markers$cell),
                            cluster_id = c('10', '11', '0', '6', '1', '2',
                                           '3', '5', '4', '9', '8', '7'))
total_markers = c()
deg_markers = c()
high_deg_markers = c()
uq_deg_markers = c()
high_uq_deg_markers = c()
deg_vs_all = c()
uq_vs_deg = c()
uq_high_vs_all_high = c()

for (i in 1:12) {
  this_cell = marker_enrich$cell[i]
  this_cluster = marker_enrich$cluster_id[i]
  markers = all.uq.markers$gene[all.uq.markers$cluster == this_cluster]
  total.n = length(markers)
  total.N = nrow(res_tpm)
  uq.deg.n = ifelse(this_cell %in% unique(filt_ann$cell),
                          table(filt_ann$cell)[this_cell],
                          0)
  uq.deg.N = nrow(filt_ann)
  uq.high.deg.n = ifelse(this_cell %in% unique(sorted_rows$celltype),
                         table(sorted_rows$celltype)[this_cell],
                         0)
  uq.high.deg.N = nrow(sorted_rows)
  deg.n = sum(degs_rnaseq$Symbol %in% markers)
  deg.N = nrow(degs_rnaseq)
  high.deg.n = sum(high_degs_rnaseq$Symbol %in% markers)
  high.deg.N = nrow(high_degs_rnaseq)
  total_markers = c(total_markers, total.n)
  deg_markers = c(deg_markers, deg.n)
  high_deg_markers = c(high_deg_markers, high.deg.n)
  uq_deg_markers = c(uq_deg_markers, uq.deg.n)
  high_uq_deg_markers = c(high_uq_deg_markers, uq.high.deg.n)
  deg_matrix <- matrix(c(deg.n, total.n - deg.n,
                         deg.N, total.N - deg.N), nrow=2, byrow=T)
  uq_vs_deg_matrix <- matrix(c(uq.deg.n, deg.n - uq.deg.n,
                             uq.deg.N, deg.N - uq.deg.N), nrow=2, byrow=T)
  uq_high_vs_deg_high_matrix <- matrix(c(uq.high.deg.n, high.deg.n - uq.high.deg.n,
                              uq.high.deg.N, high.deg.N - uq.high.deg.N), nrow=2, byrow=T)
  deg_vs_all = c(deg_vs_all, fisher.test(deg_matrix, alternative = "greater")$p.value)
  uq_vs_deg = c(uq_vs_deg, fisher.test(uq_vs_deg_matrix, alternative = "greater")$p.value)
  uq_high_vs_all_high = c(uq_high_vs_all_high, fisher.test(uq_high_vs_deg_high_matrix, alternative = "greater")$p.value)
}
marker_enrich$total_markers = total_markers
marker_enrich$deg_markers = deg_markers
marker_enrich$high_deg_markers = high_deg_markers
marker_enrich$uq_deg_markers = uq_deg_markers
marker_enrich$high_uq_deg_markers = high_uq_deg_markers
marker_enrich$deg_vs_all = deg_vs_all
marker_enrich$uq_vs_deg = uq_vs_deg
marker_enrich$uq_high_vs_all_high = uq_high_vs_all_high

marker_enrich
write.table(marker_enrich, file='cell_type_enrich_stats.tsv', sep='\t',
            row.names=F, quote=F)

marker_enrich$non_deg <- marker_enrich$total_markers - marker_enrich$deg_markers
marker_enrich$non_uq_deg <- marker_enrich$deg_markers - marker_enrich$uq_deg_markers
enrich_melt <- reshape2::melt(marker_enrich, id.vars = c('cell'),
            measure.vars=c('non_deg', 'non_uq_deg', 'uq_deg_markers'))
enrich_melt$pct = enrich_melt$value/rep(marker_enrich$total_markers, 3)
ggplot(enrich_melt, aes(x=cell, y=pct, fill=variable)) + 
  geom_bar(stat='identity', col='black') +
  theme_bw() + scale_fill_npg() + coord_flip()

circle_enrich <- reshape2::melt(marker_enrich, id.vars=c('cell'),
    measure.vars=c('deg_markers', 'uq_deg_markers', 'high_uq_deg_markers'))
circle_enrich$pct <- circle_enrich$value/c(marker_enrich$total_markers, 
            marker_enrich$deg_markers, marker_enrich$high_deg_markers) * 100
circle_enrich$pvalue = -log10(c(marker_enrich$deg_vs_all, 
                                marker_enrich$uq_vs_deg,
                                marker_enrich$uq_high_vs_all_high))

circle_enrich$significance = factor(ifelse(circle_enrich$pvalue > 5, "p < 1e-5",
                  ifelse(circle_enrich$pvalue > 3, "p < 0.001", "p > 0.001")),
                  levels = c("p > 0.001", "p < 0.001", "p < 1e-5"))

head(circle_enrich)
#circle_plot <- ggplot(circle_enrich, aes(x=cell, y=variable, size=pct)) + 
#  geom_point(aes(fill=pvalue), pch=21, col='black') + 
#  scale_fill_gradientn(colours=white_red_colors) +
#  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1),
#                     legend.position='bottom')


circle_plot <- ggplot(circle_enrich, aes(x=cell, y=variable, size=pct)) + 
  geom_point(aes(fill=significance), pch=21, col='black') + 
  scale_fill_manual(values=c("#fef0d9", "#fc8d59", "#d7301f")) +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1),
                     panel.grid=element_blank(),
                     plot.margin = margin(0.6, 2, 0.6, 0.2, "cm"))

print(circle_plot)

circle_plot_count <- ggplot(circle_enrich, aes(x=cell, y=variable, size=value)) + 
  geom_point(aes(fill=pvalue), pch=21, col='black') + 
  scale_fill_gradientn(colours=white_red_colors) +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))

# Enrichment complex heatmap
enrich_matrix <- -log10(as.matrix(marker_enrich[, 8:9]))
rownames(enrich_matrix) = marker_enrich$cell
colnames(enrich_matrix) = colnames(marker_enrich)[8:9]

gene_count_matrix <- as.matrix(marker_enrich[, c(11, 12, 6)])
gene_count_matrix <- gene_count_matrix/rowSums(gene_count_matrix)
rownames(gene_count_matrix) <- marker_enrich$cell

white_red_colors <- colorRampPalette(c('white', 'red'))(100)

enrich_ha <- ComplexHeatmap::rowAnnotation(marker_pct = 
      ComplexHeatmap::anno_barplot(gene_count_matrix, gp = gpar(fill = 2:4), 
      bar_width = 1, width = unit(13.5, "cm")))

s4_ch <- ComplexHeatmap::Heatmap(enrich_matrix, name = "mat",
                        cluster_rows=F, cluster_columns = F,
#                        col = colorRamps::matlab.like2(100),
                        col = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd"))(100), 
                        rect_gp = gpar(col = "black", lwd = 1),
                        row_names_side = "left",
                        right_annotation = enrich_ha,
                        width=unit(2, "cm"), height=unit(9.5, "cm"),
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 10),
                        column_names_side = "top")
s4_ch

top <- plot_grid(s1, s3_ch, ncol=2, rel_widths = c(1.1, 1))
bottom <- plot_grid(as.ggplot(s4_ch), plot_spacer() + theme_minimal(), rel_widths=c(1.1, 1))

cairo_pdf('Figure_4_revamp2.pdf', width=21, height=15)
plot_grid(top, bottom, nrow=2, rel_heights = c(1, 0.6))
dev.off()

first_2col <- plot_grid(s1, as.ggplot(s4_ch), ncol=1, rel_heights=c(1, 0.8))
cairo_pdf('Figure_4_revamp3.pdf', width=21, height=13)
plot_grid(first_2col, s3_ch, nrow=1, rel_widths = c(1.1, 1))
dev.off()

second_2col <- plot_grid(s1, circle_plot, ncol=1, rel_heights=c(1, 0.36))
cairo_pdf('Figure_4_revamp4.pdf', width=21, height=11.4)
plot_grid(second_2col, s3_ch, nrow=1, rel_widths = c(1.1, 1))
dev.off()


######################### PART7. Metabolic network ########################

library(mwcsr)
library(gatom)
library(data.table)
library(igraph)

res_df        <- as.data.frame(res_full)
res_df$ID     <- rownames(res_df)
de_gatom      <- res_df[,c(7,6,2,1)]
colnames(de_gatom) <- c("ID","pval","log2FC","baseMean")
rownames(de_gatom) <- 1:nrow(de_gatom)
de_gatom  <- de_gatom[complete.cases(de_gatom),]

load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rda"))
load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Hs.eg.gatom.anno.rda"))
load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rda"))

g <- makeAtomGraph(network = network,org.gatom.anno = org.Hs.eg.gatom.anno,
                   gene.de = de_gatom,met.db = met.kegg.db,met.de = NULL)
gs1 <- scoreGraph(g, k.gene = 50, k.met = NULL)
gs2 <- scoreGraph(g, k.gene = 75, k.met = NULL)
gs3 <- scoreGraph(g, k.gene = 150, k.met = NULL)

vhsolver <- virgo_solver(cplex_dir=NULL)               ## heuristic solver, needs Java 11 (conda install openjdk=11)
vsolver  <- virgo_solver(cplex_dir="/home/jovyan/bin") ## exact solver

res1 <- solve_mwcsp(vhsolver, gs1)
res2 <- solve_mwcsp(vhsolver, gs2)

res1a <- solve_mwcsp(vsolver, gs1)
res2a <- solve_mwcsp(vsolver, gs2)

saveModuleToPdf(res1$graph, file="Skaters_k50.pdf", name="Skaters_k50", n_iter=100, force=1e-5, seed=1)
saveModuleToPdf(res2$graph, file="Skaters_k75.pdf", name="Skaters_k75", n_iter=100, force=1e-5, seed=1)

saveModuleToPdf(res1a$graph, file="Skaters_k50_opt.pdf", name="Skaters_k50", n_iter=100, force=1e-5, seed=1)
saveModuleToPdf(res2a$graph, file="Skaters_k75_opt.pdf", name="Skaters_k75", n_iter=100, force=1e-5, seed=1)

res1a.ext <- addHighlyExpressedEdges(res1a$graph, gs1)
res2a.ext <- addHighlyExpressedEdges(res2a$graph, gs2)

saveModuleToPdf(res1a.ext, file="Skaters_k50_ext.pdf", name="Skaters_k50", n_iter=1000, force=1e-5, seed=1)
saveModuleToPdf(res2a.ext, file="Skaters_k75_ext.pdf", name="Skaters_k75", n_iter=1000, force=1e-5, seed=1)

save.image(file = "skaters.RData")
## load("skaters.RData")

#### heatmap of DE transporters 

res_cell <- res_tpm
res_cell$cell <- "Unknown"
res_cell$cell[res_cell$Symbol %in% pbmc.markers[pbmc.markers$cluster == "1",]$gene] <- "CD4_naive_Tcell"
res_cell$cell[res_cell$Symbol %in% pbmc.markers[pbmc.markers$cluster == "2",]$gene] <- "CD4_memory_Tcell"
res_cell$cell[res_cell$Symbol %in% pbmc.markers[pbmc.markers$cluster == "3",]$gene] <- "CD8_Tcell"
res_cell$cell[res_cell$Symbol %in% pbmc.markers[pbmc.markers$cluster == "4",]$gene] <- "natural_killer"
res_cell$cell[res_cell$Symbol %in% pbmc.markers[pbmc.markers$cluster == "5",]$gene] <- "Bcell"
res_cell$cell[res_cell$Symbol %in% pbmc.markers[pbmc.markers$cluster == "6",]$gene] <- "CD16_monocyte"
res_cell$cell[res_cell$Symbol %in% pbmc.markers[pbmc.markers$cluster == "7",]$gene] <- "dendritic_cell"
res_cell$cell[res_cell$Symbol %in% pbmc.markers[pbmc.markers$cluster == "8",]$gene] <- "platelet"
res_cell$cell[res_cell$Symbol %in% pbmc.markers[pbmc.markers$cluster == "0",]$gene] <- "CD14_monocyte"


res_cell$cell[res_cell$Symbol %in% wb.markers[wb.markers$cluster == "0",]$gene] <- "erythrocyte"
res_cell$cell[res_cell$Symbol %in% wb.markers[wb.markers$cluster == "6",]$gene] <- "neutrophil"
table(res_cell$cell)

slc_cell <- res_cell[grepl("^SLC",res_cell$Symbol),]
slc_cell <- slc_cell[slc_cell$padj <= 0.1 & slc_cell$Mean_TPM >= 10,]

slc_exp <- top18k_ann[top18k_ann$Symbol %in% slc_cell$Symbol,]
rownames(slc_exp) <- slc_exp$Symbol
slc_exp[,1:3] <- NULL
slc_exp <- slc_exp[,rownames(sorted_cond)]

slc_rows <- slc_cell[,c(2,12,5,9)]
rownames(slc_rows) <- slc_rows$Symbol
slc_rows$Symbol <- NULL

## Order by cell type and then by log fold change
slc_rows <- slc_rows[order(slc_rows$cell,slc_rows$log2FoldChange,decreasing = T),]
slc_exp <- slc_exp[rownames(slc_exp),]
write.table(slc_cell,"SLC_celltype_de.tsv",quote = F,row.names = F,sep = "\t")

## Adjust dataframes a bit more 
slc_tpms <- slc_rows[,3,drop = F]
slc_tpms$Symbol <- rownames(slc_tpms)
slc_tpms[,2:3] <- NULL
colnames(slc_tpms) <- "celltype"

sorted_tpms$Exp <- "TPM > 10"
sorted_tpms$Exp[sorted_tpms$Mean_TPM >= 50] <- "TPM > 50"
sorted_tpms$Exp[sorted_tpms$Mean_TPM >= 100] <- "TPM > 100"
sorted_tpms$Exp[sorted_tpms$Mean_TPM >= 500] <- "TPM > 500"

table(sorted_tpms$Exp)
sorted_tpms$Exp <- factor(sorted_tpms$Exp,levels = c("TPM > 10","TPM > 50","TPM > 100","TPM > 500"))
sorted_tpms$Symbol <- factor(sorted_tpms$Symbol,levels = rev(sorted_tpms$Symbol)) 



as.ggplot(pheatmap(marker_exp,cluster_rows = F,annotation_col = sorted_cond,annotation_row = sorted_rows, 
                         cluster_cols = F,scale = "row",labels_col = "", fontsize_row = 9,
                         legend = F,annotation_legend = T,border_color = "white", gaps_row = c(27,34,50,56,57,58,59,62,71,72),
                         annotation_names_row = F, annotation_names_col = F,
                         annotation_colors = list(exercise = c(after = "#FFD800",before = "#587058"))))
