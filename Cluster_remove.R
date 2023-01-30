# Author: Arboit Lorenzo
#Date: August, 24th 2021

# this CODE has to be run after
# all.clean.2 <- subset(all, subset = nFeature_RNA > 200 & nFeature_RNA <  Gene_upper[[1]] & percent.mt < 10)
# all.clean.2 <- NormalizeData(all.clean.2)
# all.clean.2 <- FindVariableFeatures(all.clean.2)
# all.clean.2 <- RunFastMNN(object.list = SplitObject(all.clean.2, split.by = "sampleInfo"))
# all.clean.2 <- RunUMAP(all.clean.2, reduction = "mnn", dims = 1:30)
# all.clean.2 <- FindNeighbors(all.clean.2, reduction = "mnn", dims = 1:30)
# all.clean.2 <- FindClusters(all.clean.2, verbose = FALSE,resolution = 0.1)

# To subset and remove single cluster and keep the remaining clusters for new analysis
sub_cluster <- subset(all.clean.2, idents = c(4, 5, 7), invert = TRUE)
DimPlot(sub_cluster,label =T,cols = cbPalette)+labs(title = 'Seurat_clusters FastMnn' )

# idents = number of clusters we want to get rid of
# re-run all the analysis with sub_cluster as obj

# sub_cluster <- subset(sub_cluster, subset = nFeature_RNA > 200 & nFeature_RNA <  Gene_upper[[1]] & percent.mt < 10)
# sub_cluster <- NormalizeData(sub_cluster)
# sub_cluster <- FindVariableFeatures(sub_cluster)
# sub_cluster <- RunFastMNN(object.list = SplitObject(sub_cluster, split.by = "sampleInfo"))
# sub_cluster <- RunUMAP(sub_cluster, reduction = "mnn", dims = 1:30)
# sub_cluster <- FindNeighbors(sub_cluster, reduction = "mnn", dims = 1:30)
# sub_cluster <- FindClusters(sub_cluster, verbose = FALSE,resolution = 0.1)
# DimPlot(sub_cluster,label =T,cols = cbPalette)+labs(title = 'Seurat_sub_clusters FastMnn' )
