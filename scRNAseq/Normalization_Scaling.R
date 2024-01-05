# Use created SeuratObject SeuratObject_Filtered.rds
# Normalize the feature expression measurements for each cell by the total expression, multiplies this by the scale factor 10'000, and log-transforms the result.
# Identify 2'000 genes (features) per dataset
# Shift the expression of each gene, so that the mean expression across cells is 0 + Scale the expression of each gene, so that the variance across cells is 1 + Regress out number of counts and mito percentage
# Save as SeuratObject_Scaled.rds

library(Seurat)
library(ggplot2)
library(Matrix)
seu <- readRDS("dir/SeuratObject_Filtered.rds")

# Normalization
seu <- Seurat::NormalizeData(seu,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)

# Variable features
seu <- Seurat::FindVariableFeatures(seu,
                                    selection.method = "vst",
                                    nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(seu), 10)
print(top10)

pdf("dir/VariableFeature.pdf")
vf_plot <- Seurat::VariableFeaturePlot(seu)
Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE)
dev.off()

# Scaling
seu <- Seurat::ScaleData(seu, vars.to.regress = c("percent.mito", "nCount_RNA"),
                         features = rownames(seu)))

saveRDS(seu, "dir/SeuratObject_Scaled.rds")