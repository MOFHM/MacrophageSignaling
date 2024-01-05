library(Seurat)

counts <- Read10X_h5("dir/file.h5", use.names = TRUE, unique.features = TRUE)
seu <- Seurat::CreateSeuratObject(counts = counts,
                                  project = "MBM")

saveRDS(seu, "dir/SeuratObject.rds")