# Use SeuratObject_Scaled.rds
# PCA on 2000 most varibale genes
# Elbow plot on 40 Pcs'
# UMAP on x PCs'
# Clusterin x PCs'
# Resoultion x
# FeaturePlots of markes
# Save as SeuratObject_Clustering.rds

library(Seurat)
library(ggplot2)
library(clustree)

seu <- readRDS("dir/SeuratObject_Scaled.rds")

# Run PCA
# Will only be run on the variable features (2'000)
seu <- Seurat::RunPCA(seu)

# View PCA
Seurat::DimPlot(seu, reduction = "pca")
ggsave("PCA_raw.pdf", path = "dir")

# ElbowPlot
pdf("dir/Elbow.pdf") 
plot(Seurat::ElbowPlot(seu, ndims = 40))
dev.off() 

# Run UMAP on x PCA
seu <- Seurat::RunUMAP(seu, dims = 1:x)

# View UMAP
Seurat::DimPlot(seu, reduction = "umap")
ggsave("UMAP_raw.pdf", path = "dir")

## Clustering
# FindNeghbours
# Constructs a SNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based 
# on the shared overlap in their local neighborhoods (Jaccard similarity)
seu <- Seurat::FindNeighbors(seu, dims = 1:x)

# FindClusters
# Modularity optimization techniques such as the Louvain algorithm (default) to iteratively group cells together, with the goal 
# of optimizing the standard modularity function. It contains a resolution parameter that sets the ‘granularity’ of the 
# downstream clustering, with increased values leading to a greater number of clusters.
seu <- Seurat::FindClusters(seu, resolution = seq(0.1, 1.5, by=0.1))
head(seu@meta.data)

# View how clusters sub-divide at increasing resolution
clustree::clustree(seu@meta.data[,grep("RNA_snn_res", colnames(seu@meta.data))],
                   prefix = "RNA_snn_res.")
ggsave("UMAP_clustertress.pdf", path = "dir") 

# View the UMAP coloring each cell according to a cluster id like this
Seurat::DimPlot(seu, group.by = "RNA_snn_res.x")
ggsave("UMAP_cluster_rx.pdf", path = "dir") 
seu <- Seurat::SetIdent(seu, value = seu$RNA_snn_res.x)

# Gene Expression of GeneY
Seurat::FeaturePlot(seu, 'GeneY', ncol=1)
ggsave("UMAP_cluster_GeneY.pdf", path = "dir") 
Seurat::VlnPlot(seu,
                features = 'GeneY',
                ncol = 1)
ggsave("Violin_cluster_GeneY.pdf", path = "dir") 

saveRDS(seu, "dir/SeuratObject_Clustering.rds")






