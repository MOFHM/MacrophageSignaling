# Macrophage Annotation based on CD163 expression

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

# Values
seuCD163_val <- FetchData(seu, vars = "CD163")$CD163

# Thresholds
q_seu <- quantile(seuCD163_val, probs = c(.1, .2, .3, .4, .5))
q0_seu <- 0
q1_seu <- q_seu[1]
q2_seu <- q_seu[2]
q3_seu <- q_seu[3]
q4_seu <- q_seu[4]
q5_seu <- q_seu[5]

# Histogram
binlist <- seq(0, 5, by=0.5)
pdf("CD163_Hist.pdf")
hist(seuCD163_val, prob = TRUE, col = 'peachpuff', border = 'black', xlab = 'Normalized Expression', main = 'CD163 Expression',
    breaks = binlist,
    ylim = c(0,1))
lines(density(seuCD163_val), lwd = 2, col = "chocolate3")
dev.off()

# Lists of CD163 values
seuCD163_values <- FetchData(seu, vars = "CD163")$CD163

# Empty lists for classifications
seuCD163_0 <- c()
seuCD163_1 <- c()
seuCD163_2 <- c()
seuCD163_3 <- c()
seuCD163_4 <- c()
seuCD163_5 <- c()

###########################################################################################
# Add CD163 Classifications @CD163_0 for positive cells
for (i in 1:length(seuCD163_values)) {
    if (seuCD163_values[i] > q0_seu){
        seuCD163_0 <- append(seuCD163_0, 'M2')
    } else {
        seuCD163_0 <- append(seuCD163_0, 'M1')}
}
seu$CD163_0 <- seuCD163_0

# Add CD163 Classifications @CD163_1 for positive cells with 10th percentile threshold
for (i in 1:length(seuCD163_values)) {
    if (seuCD163_values[i] > q1_seu){
        seuCD163_1 <- append(seuCD163_1, 'M2')
    } else {
        seuCD163_1 <- append(seuCD163_1, 'M1')}
}
seu$CD163_1 <- seuCD163_1

# Add CD163 Classifications @CD163_2 for positive cells with 20th percentile threshold
for (i in 1:length(seuCD163_values)) {
    if (seuCD163_values[i] > q2_seu){
        seuCD163_2 <- append(seuCD163_2, 'M2')
    } else {
        seuCD163_2 <- append(seuCD163_2, 'M1')}
}
seu$CD163_2 <- seuCD163_2

# Add CD163 Classifications @CD163_3 for positive cells with 30th percentile threshold
for (i in 1:length(seuCD163_values)) {
    if (seuCD163_values[i] > q3_seu){
        seuCD163_3 <- append(seuCD163_3, 'M2')
    } else {
        seuCD163_3 <- append(seuCD163_3, 'M1')}
}
seu$CD163_3 <- seuCD163_3

# Add CD163 Classifications @CD163_4 for positive cells with 40th percentile threshold
for (i in 1:length(seuCD163_values)) {
    if (seuCD163_values[i] > q4_seu){
        seuCD163_4 <- append(seuCD163_4, 'M2')
    } else {
        seuCD163_4 <- append(seuCD163_4, 'M1')}
}
seu$CD163_4 <- seuCD163_4

# Add CD163 Classifications @CD163_5 for positive cells with 50th percentile threshold
for (i in 1:length(seuCD163_values)) {
    if (seuCD163_values[i] > q5_seu){
        seuCD163_5 <- append(seuCD163_5, 'M2')
    } else {
        seuCD163_5 <- append(seuCD163_5, 'M1')}
}
seu$CD163_5 <- seuCD163_5
###########################################################################################

# Lists of Annotations
list_seu_CD163_0 <- setNames(seu$CD163_0, NULL)
list_seu_CD163_1 <- setNames(seu$CD163_1, NULL)
list_seu_CD163_2 <- setNames(seu$CD163_2, NULL)
list_seu_CD163_3 <- setNames(seu$CD163_3, NULL)
list_seu_CD163_4 <- setNames(seu$CD163_4, NULL)
list_seu_CD163_5 <- setNames(seu$CD163_5, NULL)

# DataFrames with Annotations
df1 <- data.frame(list_seu_CD163_0,list_seu_CD163_1,list_seu_CD163_2,list_seu_CD163_3,list_seu_CD163_4,list_seu_CD163_5)

write.csv(df1, "dir/CD163_Thresholds_Annotations.csv", row.names=TRUE)

# Number of annotated cells
table(seu$CD163_0)
table(seu$CD163_1)
table(seu$CD163_2)
table(seu$CD163_3)
table(seu$CD163_4)
table(seu$CD163_5)

# Load Hallmark collection from MSigDB
gmt <- msigdbr::msigdbr(species = "human", category = "H")

##########################
# Differential expression
##########################
# Q00
##########################
seu <- Seurat::SetIdent(seu, value = seu$CD163_0)
DF_seu_CD163_0 <- Seurat::FindAllMarkers(subset(seu, idents = c('M1', 'M2')),  min.pct = 0.1, test.use = "t",
                                   only.pos = TRUE)
DF_seu_CD163_0 <- subset(DF_seu_CD163_0, DF_seu_CD163_0$p_val_adj < 0.05)
DF_seu_CD163_0 <- subset(DF_seu_CD163_0, DF_seu_CD163_0$avg_log2FC > 0.5)
all_de_genesDF_seu_CD163_0 <- DF_seu_CD163_0$gene
'CD163' %in% all_de_genesDF_seu_CD163_0
DF_seu_CD163_0
##########################
# Q10
##########################
seu <- Seurat::SetIdent(seu, value = seu$CD163_1)
DF_seu_CD163_1 <- Seurat::FindAllMarkers(subset(seu, idents = c('M1', 'M2')),  min.pct = 0.1, test.use = "t",
                                   only.pos = TRUE)
DF_seu_CD163_1 <- subset(DF_seu_CD163_1, DF_seu_CD163_1$p_val_adj < 0.05)
DF_seu_CD163_1 <- subset(DF_seu_CD163_1, DF_seu_CD163_1$avg_log2FC > 0.5)
all_de_genesDF_seu_CD163_1 <- DF_seu_CD163_1$gene
'CD163' %in% all_de_genesDF_seu_CD163_1
DF_seu_CD163_1
##########################
# Q20
##########################
seu <- Seurat::SetIdent(seu, value = seu$CD163_2)
DF_seu_CD163_2 <- Seurat::FindAllMarkers(subset(seu, idents = c('M1', 'M2')),  min.pct = 0.1, test.use = "t",
                                   only.pos = TRUE)
DF_seu_CD163_2 <- subset(DF_seu_CD163_2, DF_seu_CD163_2$p_val_adj < 0.05)
DF_seu_CD163_2 <- subset(DF_seu_CD163_2, DF_seu_CD163_2$avg_log2FC > 0.5)
all_de_genesDF_seu_CD163_2 <- DF_seu_CD163_2$gene
'CD163' %in% all_de_genesDF_seu_CD163_2
DF_seu_CD163_2
##########################
# Q30
##########################
seu <- Seurat::SetIdent(seu, value = seu$CD163_3)
DF_seu_CD163_3 <- Seurat::FindAllMarkers(subset(seu, idents = c('M1', 'M2')),  min.pct = 0.1, test.use = "t",
                                   only.pos = TRUE)
DF_seu_CD163_3 <- subset(DF_seu_CD163_3, DF_seu_CD163_3$p_val_adj < 0.05)
DF_seu_CD163_3 <- subset(DF_seu_CD163_3, DF_seu_CD163_3$avg_log2FC > 0.5)
all_de_genesDF_seu_CD163_3 <- DF_seu_CD163_3$gene
'CD163' %in% all_de_genesDF_seu_CD163_3
DF_seu_CD163_3
##########################
# Q40
##########################
seu <- Seurat::SetIdent(seu, value = seu$CD163_4)
DF_seu_CD163_4 <- Seurat::FindAllMarkers(subset(seu, idents = c('M1', 'M2')),  min.pct = 0.1, test.use = "t",
                                   only.pos = TRUE)
DF_seu_CD163_4 <- subset(DF_seu_CD163_4, DF_seu_CD163_4$p_val_adj < 0.05)
DF_seu_CD163_4 <- subset(DF_seu_CD163_4, DF_seu_CD163_4$avg_log2FC > 0.5)
all_de_genesDF_seu_CD163_4 <- DF_seu_CD163_4$gene
'CD163' %in% all_de_genesDF_seu_CD163_4
DF_seu_CD163_4
##########################
# Q50
##########################
seu <- Seurat::SetIdent(seu, value = seu$CD163_5)
DF_seu_CD163_5 <- Seurat::FindAllMarkers(subset(seu, idents = c('M1', 'M2')),  min.pct = 0.1, test.use = "t",
                                   only.pos = TRUE)
DF_seu_CD163_5 <- subset(DF_seu_CD163_5, DF_seu_CD163_5$p_val_adj < 0.05)
DF_seu_CD163_5 <- subset(DF_seu_CD163_5, DF_seu_CD163_5$avg_log2FC > 0.5)
all_de_genesDF_seu_CD163_5 <- DF_seu_CD163_5$gene
'CD163' %in% all_de_genesDF_seu_CD163_5
DF_seu_CD163_5

##########################
# Enrichment with CD163
##########################
# Q00
##########################
# Subsets
DF_seu_CD163_M1_0 <- subset(DF_seu_CD163_0, cluster=='M1')[,7]
DF_seu_CD163_M2_0 <- subset(DF_seu_CD163_0, cluster=='M2')[,7]
M1_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_0, cluster=='M1')$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_0, cluster=='M2')$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]

##########################
# Q10
##########################
# Subsets
DF_seu_CD163_M1_1 <- subset(DF_seu_CD163_1, cluster=='M1')[,7]
DF_seu_CD163_M2_1 <- subset(DF_seu_CD163_1, cluster=='M2')[,7]
M1_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_1, cluster=='M1')$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_1, cluster=='M2')$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]

##########################
# Q20
##########################
# Subsets
DF_seu_CD163_M1_2 <- subset(DF_seu_CD163_2, cluster=='M1')[,7]
DF_seu_CD163_M2_2 <- subset(DF_seu_CD163_2, cluster=='M2')[,7]
M1_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_2, cluster=='M1')$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_2, cluster=='M2')$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]

##########################
# Q30
##########################
# Subsets
DF_seu_CD163_M1_3 <- subset(DF_seu_CD163_3, cluster=='M1')[,7]
DF_seu_CD163_M2_3 <- subset(DF_seu_CD163_3, cluster=='M2')[,7]
M1_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_3, cluster=='M1')$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_3, cluster=='M2')$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]

##########################
# Q40
##########################
# Subsets
DF_seu_CD163_M1_4 <- subset(DF_seu_CD163_4, cluster=='M1')[,7]
DF_seu_CD163_M2_4 <- subset(DF_seu_CD163_4, cluster=='M2')[,7]
M1_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_4, cluster=='M1')$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_4, cluster=='M2')$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]

##########################
# Q50
##########################
# Subsets
DF_seu_CD163_M1_5 <- subset(DF_seu_CD163_5, cluster=='M1')[,7]
DF_seu_CD163_M2_5 <- subset(DF_seu_CD163_5, cluster=='M2')[,7]
M1_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_5, cluster=='M1')$gene,
                                                universe = rownames(seu),
                                                pAdjustMethod = "BH",
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M1_enrich@result[which(M1_enrich@result$p.adjust<0.05),]

M2_enrich <- clusterProfiler::enricher(gene = subset(DF_seu_CD163_5, cluster=='M2')$gene,
                                       universe = rownames(seu),
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 0.05,
                                       TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
M2_enrich@result[which(M2_enrich@result$p.adjust<0.05),]
