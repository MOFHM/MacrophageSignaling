# Macrophage Annotation with Extended Set Markers
# Azizi et al. Cell (2018): 10.1016/j.cell.2018.05.060

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

seu <- readRDS("dir/SeuratObject_Celldex.rds")
DefaultAssay(seu) <- "RNA"
seu <- Seurat::SetIdent(seu, value = seu$orig.ident)

##########################
# Input Azizi markers
aziziM1 <- c('CCL5','CCR10','CCR7','CD40','CD74','CD80','CD86','CXCL10','CXCL11','CXCL9','FCGR1A','FCGR1B','FCGR1C','HLA-DMA','HLA-DMB','HLA-DRA','HLA-DRB1','HLA-DRB3','IDO1','IL12A','IL12B','IL1A','IL1B','IL23A','IL6','IRF1','IRF5','KYNU','NOS2','TNF')
aziziM2 <- c('ARG1','ARG2','CCL13','CCL17','CCL18','CCL20','CCL22','CCL24','CCL4','CD163','CD200R1','CD274','CD276','CLEC7A','CSF1R','CTSC','CTSA','CTSB','CTSD','EGF','FASLG','FCER2','FCGR2A','FCGR2B','FCGR2C','FN1','IL10','IL1R2','IL1RN','IL4R','IRF4','LYVE1','MARCO','MMP14','MMP19','MMP9','MRC1',
        'MSR1','PDCD1LG2','TGFB1','TGFB2','TGFB3','TNFSF12','TNFSF8','VEGFA','VEGFB','VEGFC','FIGF','VTCN1','WNT7B')
##########################

##########################
# Retain markers present in SeuratObject
seuaziziM1 <- aziziM1[aziziM1 %in% rownames(seu)]
seuaziziM2 <- aziziM2[aziziM2 %in% rownames(seu)]

##########################
# Azizi Classification @Azizi
th <- 0.0

seu <- Seurat::AddModuleScore(seu,
                              features = list(seuaziziM1),
                              name = "M1_Azizi")
seu <- Seurat::AddModuleScore(seu,
                              features = list(seuaziziM2),
                              name = "M2_Azizi")
eisseu <- seu$M1_Azizi1
zweiseu <- seu$M2_Azizi1
newseu <- c()
for (i in 1:length(colnames(seu))) {
  n1 <- eisseu[i]
  n2 <- zweiseu[i]
  if (n1 > n2){if (n1 > th){newseu <- append(newseu, 'M1')} else {newseu <- append(newseu, 'Na')}}
  if (n2 > n1){if (n2 > th){newseu <- append(newseu, 'M2')} else {newseu <- append(newseu, 'Na')}}
}
seu$Azizi <- newseu

##########################
table(seu$Azizi)

##########################
# Load Hallmark collection from MSigDB
gmt <- msigdbr::msigdbr(species = "human", category = "H")

##########################
# Azizi
# Differential expression
# and Enrichment
##########################

seu <- Seurat::SetIdent(seu, value = seu$Azizi)
DF_seu_azizi <- Seurat::FindAllMarkers(subset(seu, idents = c('M1', 'M2')),  min.pct = 0.1, test.use = "t",
                                   only.pos = TRUE)
DF_seu_azizi <- subset(DF_seu_azizi, DF_seu_azizi$p_val_adj < 0.05)
DF_seu_azizi <- subset(DF_seu_azizi, DF_seu_azizi$avg_log2FC > 0.75)
all_de_genesDF_seu_azizi <- DF_seu_azizi$gene
intersect(DF_seu_azizi$gene, c(aziziM1, aziziM2))
# Subsets
DF_seu_azizi_M1 <- noquote(subset(DF_seu_azizi, cluster=='M1')[,7])
DF_seu_azizi_M2 <- noquote(subset(DF_seu_azizi, cluster=='M2')[,7])
# Enrichment all Genes
M1_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_azizi, cluster=='M1')$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_azizi, cluster=='M2')$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]

# Enrichment new Genes
M1_enrich <- clusterProfiler::enricher(gene = subset(subset(DF_seu_azizi, cluster=='M1'), !(gene %in% aziziM1))$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(subset(DF_seu_azizi, cluster=='M2'), !(gene %in% aziziM2))$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]
