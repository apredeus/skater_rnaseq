stringsAsFactors = F
setwd("~/tmp")

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

library(gatom)
library(data.table)
library(igraph)

rsem_exp     <- read.table("fixed_rsem.genes.counts.tsv", header=T, row.names = 1, check.names=F)
dim(rsem_exp)
head(rsem_exp)
## PAR genes are removed, thus no collapsing is needed

cond          <- read.table("Conditions.txt", row.names = 1, header = T, stringsAsFactors = T)
ann           <- rsem_exp[,c(1:2)]
exp           <- rsem_exp[,c(3:16)]
colSums(exp[, c(1:14)])

dds_full      <- DESeqDataSetFromMatrix(countData=round(exp,0),colData = cond,design = ~ Donor + Cond)
dds_full$Cond <- relevel(dds_full$Cond, ref = "before")
dds_full      <- DESeq(dds_full)
res_full      <- results(dds_full, contrast=c("Cond","after","before"))
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

vhsolver <- virgo_solver(cplex_dir=NULL)               ## heuristic solver
vsolver  <- virgo_solver(cplex_dir="/home/jovyan/bin") ## exact solver

res1 <- solve_mwcsp(vhsolver, gs1)
res2 <- solve_mwcsp(vhsolver, gs2)

res1a <- solve_mwcsp(vsolver, gs1)
res2a <- solve_mwcsp(vsolver, gs2)

saveModuleToPdf(res1$graph, file="Skaters_k50.pdf", name="Skaters_k50", n_iter=100, force=1e-5, seed=1)
saveModuleToPdf(res2$graph, file="Skaters_k75.pdf", name="Skaters_k75", n_iter=100, force=1e-5, seed=1)

saveModuleToPdf(res1a$graph, file="Skaters_k50_opt.pdf", name="Skaters_k50", n_iter=100, force=1e-5, seed=1)
saveModuleToPdf(res2a$graph, file="Skaters_k75_opt.pdf", name="Skaters_k75", n_iter=100, force=1e-5, seed=1)




m <- res$graph


















solve <- sgmwcs.solver("/home/jovyan/bin/sgmwcs", nthreads = 4, timeLimit = 60) 
m <- solve(gs)

m <- solveSgmwcsRandHeur(gs, max.iterations = 2000)
gs1 <- scoreGraph(g, k.gene = 50, k.met = NULL)
m1  <- solve(gs1)
gs2 <- scoreGraph(g, k.gene = 50, k.met = NULL)
m2  <- solve2(gs2)
m3  <- solve2(gs3)
m4  <- solve2(gs4)

saveModuleToDot(m1, file="skaters.k50.v1.dot", name="skaters.k50.v1")
saveModuleToDot(m2, file="skaters.k50.v2.dot", name="skaters.k50.v2")
saveModuleToDot(m3, file="skaters.k75.v1.dot", name="skaters.k75.v1")
saveModuleToDot(m4, file="skaters.k150.v1.dot", name="skaters.k150.v1")

system("neato -Tsvg skaters.k50.v1.dot > skaters.k50.v1.svg", ignore.stderr = T)
system("neato -Tsvg skaters.k50.v2.dot > skaters.k50.v2.svg", ignore.stderr = T)
system("neato -Tsvg skaters.k75.v1.dot > skaters.k75.v1.svg", ignore.stderr = T)
system("neato -Tsvg skaters.k150.v1.dot > skaters.k150.v1.svg", ignore.stderr = T)


saveModuleToPdf(m2, file="skaters.k50.v2.pdf", name="skaters.k50.v2", n_iter=100, force=1e-5, seed=1)

m2.ext <- addHighlyExpressedEdges(m2, gs2)
saveModuleToDot(m2.ext, file="skaters.k50.v3.dot", name="skaters.k50.v3")
system("neato -Tsvg skaters.k50.v3.dot > skaters.k50.v3.svg", ignore.stderr = T)

m2.ext2 <- addHighlyExpressedEdges(m2, gs2, top=6000)
saveModuleToDot(m2.ext2, file="skaters.k50.v4.dot", name="skaters.k50.v4")
system("neato -Tsvg skaters.k50.v4.dot > skaters.k50.v4.svg", ignore.stderr = T)


