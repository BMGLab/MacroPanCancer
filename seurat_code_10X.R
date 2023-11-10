library(Seurat)
library(dplyr)
library(hdf5r)


umi_matrix.data <- Read10X(data.dir = "path/to/file")
#umi_matrix.data <- Read10X_h5("../Desktop/LAB/DATASETS/O1/GSM3319032_sample_1-1_filtered_gene_bc_matrices_h5.h5")

seurat_obj <- CreateSeuratObject(counts = umi_matrix.data, min.cells = 3, min.features = 200)                      

# Normalization
seurat_obj <- NormalizeData(seurat_obj)
dim(seurat_obj)

# Find variable features 
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Özellik sıralamasını ve boyutlandırmayı yap
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.8))
seurat_obj <- RunUMAP(seurat_obj, dims = (1:10))
DimPlot(seurat_obj, reduction = "umap")

# macrophage control
VlnPlot(seurat_obj, features = c("LYZ", "CD68"))
VlnPlot(seurat_obj, features = c("CD163", "IL4I1"))
VlnPlot(seurat_obj, features = c("CD14", "CD64"))

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster2.markers <- FindMarkers(seurat_obj, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
cluster3.markers <- FindMarkers(seurat_obj, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)
cluster10.markers <- FindMarkers(seurat_obj, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 5)

# rename clusters according to markers
new.cluster.ids <- c("0","1", "2", "3","4","5","6","7","8", "9", "Macrophages","11","12","13","14","15","16","17","18","19","20","21","22","23", "24", "25")
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

#subset
Idents(object = seurat_obj) <- seurat_obj$seurat_clusters
seurat_obj <- RenameIdents(seurat_obj, '10' = "Macrophages")
seurat_obj <- subset(x = seurat_obj, idents = "Macrophages")

# Identification of M1/M2 macrophages
## Markers for M1 macrophages
FeaturePlot(GSE130148_Lung, features = c("CD80", "CD86"), label = F)
## Markers for M2 macrophages
FeaturePlot(GSE130148_Lung, features = c("CD163", "CD206", "CD68"), label = F)

VlnPlot(seurat_obj, features = c("LYZ", "IL4I1"))
VlnPlot(seurat_obj, features = c("CD163", "CD206", "CD68"))
VlnPlot(seurat_obj, features = c("CD14", "CD16", "CD64"))

dim(seurat_obj)
saveRDS(seurat_obj, file = "P1_seurat.RDS")
