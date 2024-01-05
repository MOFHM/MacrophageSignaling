# Macrophage Annotation with Proteomics Markers

library(Seurat)
library(edgeR)
library(limma)
library(dplyr)
library(dittoSeq)
library(clustree)
library(celldex)
library(data.table)
library(rstatix)
library(org.Hs.eg.db)
library(ggnewscale)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggpubr)

seu <- readRDS("dir/SeuratObject__Celldex.rds")
DefaultAssay(seu) <- "RNA"
seu <- Seurat::SetIdent(seu, value = seu$orig.ident)

##########################
# Input proteomics markers
protm1 <- read.csv("Proteomics_Signature M1.csv",header=F)$V1
protm2 <- read.csv("Proteomics_Signature M2.csv",header=F)$V1

##########################
# Retain markers present in SeuratObject
seuprotm1 <- protm1[protm1 %in% rownames(seu)]
seuprotm2 <- protm2[protm2 %in% rownames(seu)]

##########################
# Proteomics Classification @Proteomics
th <- 0.0              

seu <- Seurat::AddModuleScore(seu,
                              features = list(seuprotm1),
                              name = "M1_Proteomics")
seu <- Seurat::AddModuleScore(seu,
                              features = list(seuprotm2),
                              name = "M2_Proteomics")
eisseu <- seu$M1_Proteomics1
zweiseu <- seu$M2_Proteomics1

newseu <- c()
for (i in 1:length(colnames(seu))) {
  n1 <- eisseu[i]
  n2 <- zweiseu[i]
  if (n1 > n2){if (n1 > th){newseu <- append(newseu, 'M1')} else {newseu <- append(newseu, 'Na')}}
  if (n2 > n1){if (n2 > th){newseu <- append(newseu, 'M2')} else {newseu <- append(newseu, 'Na')}}
}
seu$Proteomics <- newseu

table(seu$Proteomics)

# Check presence
# Split
seu_list <- Seurat::SplitObject(seu, split.by = "orig.ident")
for (i in 1:length(seu_list)) {
  print('-----')
  print(table(seu_list[[i]]@meta.data$orig.ident))
  print(table(seu_list[[i]]@meta.data$Proteomics))
}

# Load Hallmark collection from MSigDB
gmt <- msigdbr::msigdbr(species = "human", category = "H")

##########################
# Proteomics
# Differential expression
# and Enrichment
##########################

seu <- Seurat::SetIdent(seu, value = seu$Proteomics)
DF_seu_prot <- Seurat::FindAllMarkers(subset(seu, idents = c('M1', 'M2')),  min.pct = 0.1, test.use = "t",
                                   only.pos = TRUE)
DF_seu_prot <- subset(DF_seu_prot, DF_seu_prot$p_val_adj < 0.05)
DF_seu_prot <- subset(DF_seu_prot, DF_seu_prot$avg_log2FC > 0.75)
all_de_genesDF_seu_prot <- DF_seu_prot$gene
intersect(DF_seu_prot$gene, c(protm1, protm2))
# Subsets
DF_seu_prot_M1 <- noquote(subset(DF_seu_prot, cluster=='M1')[,7])
DF_seu_prot_M2 <- noquote(subset(DF_seu_prot, cluster=='M2')[,7])

# Enrichment all Genes
M1_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_prot, cluster=='M1')$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_prot, cluster=='M2')$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]

# Enrichment new Genes
M1_enrich <- clusterProfiler::enricher(gene = subset(subset(DF_seu_prot, cluster=='M1'), !(gene %in% protm1))$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(subset(DF_seu_prot, cluster=='M2'), !(gene %in% protm2))$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]
