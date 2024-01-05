# Read the previously with python converted .h5ad files (already converted to .tsv.gz files) and create a SeuratObject out of it

library(Seurat)
library(Matrix)

sparse_matrix <- Seurat::Read10X(data.dir = dir)
metadata <- read.csv('metadata.csv')    
rownames(metadata) <- metadata$index

seu <- Seurat::CreateSeuratObject(counts = sparse_matrix, meta.data = metadata)                                 
saveRDS(seu, "dir/SeuratObject.rds")
