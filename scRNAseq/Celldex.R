# Use SeuratObject_Clustering.rds
# Exclude clusters X, Y, Z (tumor cells)
# Reanotate using celldex and SingleR
# Save as SeuratObject_Celldex.rds

library(Seurat)
library(ggplot2)
library(celldex)
library(SingleR)

## Data input
seu <- readRDS("dir/SeuratObject_Clustering.rds")

# Check presence of certain datasets in clusters
# Split
seu_list <- Seurat::SplitObject(seu, split.by = "RNA_snn_res.x")
for (i in 1:length(seu_list)) {
  print('-----')
  # Dataset
  print(table(seu_list[[i]]@meta.data$patientID))
  # Cluster
  print(table(Idents(seu_list[[i]])))
}

# Specify clusters of stromal cells to include
incl <- c(A, B, C)
# Split into subsets
seu <- subset(seu, idents = incl)

# Seu organization
DefaultAssay(seu) <- "RNA"
Idents(seu) <- "patientID"
print(seu)

ref <- BlueprintEncodeData()
class(ref)
table(ref$label.main)

# Define which cell types to include
ye <- c('B-cells', 'CD4+ T-cells', 'CD8+ T-cells', 'DC', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Macrophages', 'Monocytes', 'NK cells')  
ref <- ref[,ref$label.main %in% ye]
table(ref$label.main)

seu_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu, slot = "data"),
                                    ref = ref,
                                    labels = ref$label.main)
head(seu_SingleR)

SingleR::plotDeltaDistribution(seu_SingleR, dots.on.top = FALSE)
ggsave("SingleR_Distribution.pdf", path = "dir") 

singleR_labels <- seu_SingleR$labels
seu$SingleR_annot <- singleR_labels

Seurat::DimPlot(seu, group.by = "SingleR_annot")
ggsave("DimPlot_annot.pdf", path = "dir") 

# Check presence of certain datasets in cell type classes
# Split
seu_list <- Seurat::SplitObject(seu, split.by = "SingleR_annot")
for (i in 1:length(seu_list)) {
  print('-----')
  print(table(seu_list[[i]]@meta.data$SingleR_annot))
  # Dataset
  print(table(seu_list[[i]]@meta.data$patientID))
  # Cluster
  print(table(Idents(seu_list[[i]])))
}

# Macrophages
seu <- subset(seu, SingleR_annot == 'Macrophages')
head(seu@meta.data)

# Save
saveRDS(seu, "dir/SeuratObject_Celldex.rds")