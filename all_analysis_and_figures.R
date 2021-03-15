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

ery <- wb.markers[wb.markers$cluster == "0",]
ery <- ery[! grepl("^HB",ery$gene),] ## lose the globins since we depleted them
ery$cell <- "erythrocyte"
ery <- ery[1:20,]

neutro <- wb.markers[wb.markers$cluster == "6",]
neutro$cell <- "neutrophil"
neutro <- neutro[1:20,]

### Let's process 10X PBMC10k dataset the same way and get all other markers 
download.file("https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_filtered_feature_bc_matrix.h5",
              destfile = "pbmc10k_filt.h5")
filt.matrix <- Read10X_h5("pbmc10k_filt.h5",use.names = T)
pbmc        <- CreateSeuratObject(counts = filt.matrix)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 15) ## crude filtering for doublets
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc    <- NormalizeData(pbmc)
pbmc    <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

## Again, we want courser clustering - overall B/T/Monocyte/NK/platelet 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.2)

DimPlot(pbmc,label.size = 6,repel = T,label = T)
FeaturePlot(pbmc,c("LILRA4"))
FeaturePlot(pbmc,c("PPBP")) ## got to split clusters, these are DCs and platelets together

plat_barcodes <- colnames(subset(x = pbmc, subset = PPBP > 2))
pbmc@meta.data$seurat_clusters <- as.character(pbmc@meta.data$seurat_clusters)
pbmc@meta.data$seurat_clusters[rownames(pbmc@meta.data) %in% plat_barcodes] <- '8'
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

cd14_mono <- pbmc.markers[pbmc.markers$cluster == "0",]
cd14_mono$cell <- "CD14_monocyte"
cd14_mono <- cd14_mono[1:20,]

cd4_naive <- pbmc.markers[pbmc.markers$cluster == "1",]
cd4_naive$cell <- "CD4_naive_Tcell"
cd4_naive <- cd4_naive[1:20,]

cd4_memory <- pbmc.markers[pbmc.markers$cluster == "2",]
cd4_memory$cell <- "CD4_memory_Tcell"
cd4_memory <- cd4_memory[1:20,]

cd8_tcell <- pbmc.markers[pbmc.markers$cluster == "3",]
cd8_tcell$cell <- "CD8_Tcell"
cd8_tcell <- cd8_tcell[1:20,]

nk <- pbmc.markers[pbmc.markers$cluster == "4",]
nk$cell <- "natural_killer"
nk <- nk[1:20,]

bcell <- pbmc.markers[pbmc.markers$cluster == "5",]
bcell$cell <- "Bcell"
bcell <- bcell[1:20,]

cd16_mono <- pbmc.markers[pbmc.markers$cluster == "6",]
cd16_mono$cell <- "CD16_monocyte"
cd16_mono <- cd16_mono[1:20,]

dc <- pbmc.markers[pbmc.markers$cluster == "7",]
dc$cell <- "dendritic_cell"
dc <- dc[1:20,]

plat <- pbmc.markers[pbmc.markers$cluster == "8",]
plat$cell <- "platelet"
plat <- plat[1:20,]

all_top20_markers <- rbind(ery,neutro,cd14_mono,cd16_mono,cd4_naive,cd4_memory,cd8_tcell,nk,bcell,dc,plat)
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

## Let's take expression matrix after donor correction, and find all the defined scRNAseq cell-type specific markers. 
cbat_select <- top18k_ann[top18k_ann$Symbol %in% all_top20_markers$gene,]
dim(cbat_select)
cbat_markers <- merge(all_top20_markers[,7:8],cbat_select,by.x = "gene",by.y = "Symbol")
cbat_markers$Ensembl <- NULL
cbat_markers$Gene_type <- NULL
markers_melt <- melt(cbat_markers)
markers_melt$exercise <- gsub("\\d","",markers_melt$variable,perl = T)
markers_melt$exercise <- factor(markers_melt$exercise,levels = c("before","after"))
markers_melt$variable <- NULL

## boxplot of all markers, before and after exercise
r1 <- ggplot(markers_melt,aes(x=cell,y=value,fill=exercise)) + geom_boxplot() + theme_bw() + 
  scale_fill_manual(values = c("#587058","#FFD800")) + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  ylab("Top 20 single cell markers") + xlab("cell type")

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

r2 <- wrap_elements(grid::textGrob('Unchanged')) + hmap[["erythrocyte"]] + 
  hmap[["CD16_monocyte"]] + hmap[["Bcell"]] + hmap[["CD4_memory_Tcell"]] + plot_layout(ncol = 5, widths = c(1,3,3,3,3))
r3 <- wrap_elements(grid::textGrob('Decreased')) + hmap[["CD4_naive_Tcell"]] + 
  hmap[["CD8_Tcell"]] + hmap[["natural_killer"]] + hmap[["dendritic_cell"]] + plot_layout(ncol = 5, widths = c(1,3,3,3,3))
r4 <- wrap_elements(grid::textGrob('Increased')) + hmap[["neutrophil"]] + 
  hmap[["CD14_monocyte"]] + hmap[["platelet"]] + plot_spacer() + plot_layout(ncol = 5, widths = c(1,3,3,3,3))


## Patchwork figure 3, 15 x 18. TODO: sort cells by unchanged/dec/inc. 
r1/r2/r3/r4

################################## PART6. Annotate and plot more "array-filtered" genes ########################

## This is a logical continuation of PART3, but with blackjack and single cell markers. 
dim(filtered)
head(filtered)

### We can make a volcano plot of the filtered genes. It's not super informative though. 

filt_label <- subset(filtered, Mean_TPM > 10 & Gene_type == "protein_coding" & abs(log2FoldChange) >= 0.5)
write.table(filt_label,"Array_filtered_protein_tpm10_de.tsv",quote = F,row.names = F,sep = "\t")

ggplot(filtered, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant, size = log10(Mean_TPM+1)), stroke = 0) + 
  geom_text_repel(data = filt_label, aes(label = Symbol), size = 2) + 
  scale_size_continuous(range = c(0.2,3)) + 
  scale_x_continuous(limits = c(-2, 2)) + 
  scale_color_manual(values = c("red", "gray80")) + 
  xlab("log fold change") + 
  ylab("-log10 (adjusted p-value)") + 
  ggtitle("Volcano plot, before vs. after exercise") + theme_bw() + 
  theme(panel.grid.major = element_line(size = 0.25))

table(all_top20_markers$cell)
filt_ann <- filtered
filt_ann$cell <- "Unknown"
filt_ann$cell[filt_ann$Symbol %in% pbmc.markers[pbmc.markers$cluster == "1",]$gene] <- "CD4_naive_Tcell"
filt_ann$cell[filt_ann$Symbol %in% pbmc.markers[pbmc.markers$cluster == "2",]$gene] <- "CD4_memory_Tcell"
filt_ann$cell[filt_ann$Symbol %in% pbmc.markers[pbmc.markers$cluster == "3",]$gene] <- "CD8_Tcell"
filt_ann$cell[filt_ann$Symbol %in% pbmc.markers[pbmc.markers$cluster == "4",]$gene] <- "natural_killer"
filt_ann$cell[filt_ann$Symbol %in% pbmc.markers[pbmc.markers$cluster == "5",]$gene] <- "Bcell"
filt_ann$cell[filt_ann$Symbol %in% pbmc.markers[pbmc.markers$cluster == "6",]$gene] <- "CD16_monocyte"
filt_ann$cell[filt_ann$Symbol %in% pbmc.markers[pbmc.markers$cluster == "7",]$gene] <- "dendritic_cell"
filt_ann$cell[filt_ann$Symbol %in% pbmc.markers[pbmc.markers$cluster == "8",]$gene] <- "platelet"
filt_ann$cell[filt_ann$Symbol %in% pbmc.markers[pbmc.markers$cluster == "0",]$gene] <- "CD14_monocyte"


filt_ann$cell[filt_ann$Symbol %in% wb.markers[wb.markers$cluster == "0",]$gene] <- "erythrocyte"
filt_ann$cell[filt_ann$Symbol %in% wb.markers[wb.markers$cluster == "6",]$gene] <- "neutrophil"
table(filt_ann$cell)

filt_label2 <- subset(filt_ann, Mean_TPM > 10 & Gene_type == "protein_coding" & abs(log2FoldChange) >= 0.5)
write.table(filt_label2,"Array_filtered_celltype_de.tsv",quote = F,row.names = F,sep = "\t")
