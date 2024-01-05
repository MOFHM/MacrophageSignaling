# Macrophage Annotation with Core Set Markers
# Murray et al. Immunity (2014): 10.1016/j.immuni.2014.06.008

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
# Input Murray markers
murrayM1 <- c('CCL5','CCL18','CCR7','CD40','CXCL10','CXCL11','CXCL9','GBP1','IDO1','IL12A','IL12B','IL1B','IL23A','IL6','CXCL8','IRF1','IRF5','KYNU','PTX3','TNF')
murrayM2 <- c('ADORA3','ALOX15','CCL13','CCL17','CCL18','CCL4','CD163','CD200R1','F13A1','FN1','GATA3','IL17RB','IL4R','IRF4','MARCO','MMP1','MMP12','MRC1','SOCS1','SOCS3','STAB1','TG','TGFB1','TGFBR2','TGM2') 
##########################

##########################
# Retain markers present in SeuratObject
seumurrayM1 <- murrayM1[murrayM1 %in% rownames(seu)]
seumurrayM2 <- murrayM2[murrayM2 %in% rownames(seu)]

##########################
# Murray Classification @Murray
th <- 0.0

seu <- Seurat::AddModuleScore(seu,
                              features = list(seumurrayM1),
                              name = "M1_Murray")
seu <- Seurat::AddModuleScore(seu,
                              features = list(seumurrayM2),
                              name = "M2_Murray")
eisseu <- seu$M1_Murray1
zweiseu <- seu$M2_Murray1
newseu <- c()
for (i in 1:length(colnames(seu))) {
  n1 <- eisseu[i]
  n2 <- zweiseu[i]
  if (n1 > n2){if (n1 > th){newseu <- append(newseu, 'M1')} else {newseu <- append(newseu, 'Na')}}
  if (n2 > n1){if (n2 > th){newseu <- append(newseu, 'M2')} else {newseu <- append(newseu, 'Na')}}
}
seu$Murray <- newseu

##########################
table(seu$Murray)

##########################
# Load Hallmark collection from MSigDB
gmt <- msigdbr::msigdbr(species = "human", category = "H")
##########################

##########################
# Murray
# Differential expression
# and Enrichment
##########################

seu <- Seurat::SetIdent(seu, value = seu$Murray)
DF_seu_murray <- Seurat::FindAllMarkers(subset(seu, idents = c('M1', 'M2')),  min.pct = 0.1, test.use = "t",
                                   only.pos = TRUE)
DF_seu_murray <- subset(DF_seu_murray, DF_seu_murray$p_val_adj < 0.05)
DF_seu_murray <- subset(DF_seu_murray, DF_seu_murray$avg_log2FC > 0.75)
all_de_genesDF_seu_murray <- DF_seu_murray$gene
intersect(DF_seu_murray$gene, c(murrayM1, murrayM2))
# Subsets
DF_seu_murray_M1 <- noquote(subset(DF_seu_murray, cluster=='M1')[,7])
DF_seu_murray_M2 <- noquote(subset(DF_seu_murray, cluster=='M2')[,7])
# Enrichment all Genes
M1_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_murray, cluster=='M1')$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_murray, cluster=='M2')$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]

# Enrichment new Genes
M1_enrich <- clusterProfiler::enricher(gene = subset(subset(DF_seu_murray, cluster=='M1'), !(gene %in% murrayM1))$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(subset(DF_seu_murray, cluster=='M2'), !(gene %in% murrayM2))$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]
