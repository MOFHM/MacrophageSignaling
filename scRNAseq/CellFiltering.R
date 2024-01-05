# Use created SeuratObject SeuratObject.rds
# Calculate mitochondrial genes in %
# Cell filtering (Genes in 3 cells, 500-9000 genes per cell, <10% mitochondrial genes)
# Evaluate with ViolinPlot
# Save as SeuratObject_Filtered.rds

library(Seurat)
library(ggplot2)
library(Matrix)
seu <- readRDS("dir/SeuratObject.rds")

seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^MT-", 
                                    col.name = "percent.mito")

filter_num_expressed_cells <- rownames(seu)[Matrix::rowSums(seu) > 3]
seu <- subset(seu, subset = nFeature_RNA > 500 & 
                nFeature_RNA < 9000 &
                percent.mito < 10)
seu <- subset(seu, features = filter_num_expressed_cells)
dim(seu)

pdf("dir/Violin_nFeature_Mito.pdf")
plot(Seurat::VlnPlot(seu, features = c("nFeature_RNA",
                                       "percent_mito")))
dev.off()

saveRDS(seu, "dir/SeuratObject_Filtered.rds")