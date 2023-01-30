# Author: Arboit Lorenzo
#Date: August, 24th 2021

# Load libraries
library(scater)
library(Seurat)
library(dplyr)
library(hdf5r)
library(clustree)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(monocle3)
library(magrittr)
library(patchwork)
library(ggraph)
library(SeuratWrappers)
library(ggplot2)
library(Signac)


cbPalette <- c("#000000","#009292","#ff6db6","#ffb6db","#004949",
               "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
               "#920000","#924900","#db6d00","#24ff24","#ffff6d")

targetGenes=unique(c("Insert genes of interest as their HGNC code"))

# Allow to allocate more momory -> solution to error vector allocation
invisible(utils::memory.limit(64000))

###############################
# DATA IMPORT FROM 10X GENOMICS
###############################

# Import data for CONTROL:
# How to read in 10X data for a single sample (output is a sparse matrix)
ctrol_counts=Read10X_h5('Data_name.h5', use.names = TRUE, unique.features = TRUE)

# Turn count matrix into a Seurat object (output is a Seurat object). We decide to keep only the genes that are present in at least 3 cells (min.cells)
# and the cells that express at least 200 genes (min.features), avoid cells with features <200, they could be death cells
ctrol_obj <- CreateSeuratObject(counts = ctrol_counts, project = "CTRL", min.cells = 3, min.features = 200)
ctrol_obj.totalCell=as.integer(dim(ctrol_counts)[2])

# Import data for MUTANT
exp_counts=Read10X_h5('Dara_name.h5', use.names = TRUE, unique.features = TRUE)
exp_obj.totalCell=as.integer(dim(exp_counts)[2])
exp_obj <- CreateSeuratObject(counts = exp_counts, project = "MUT", min.cells = 3, min.features = 200)

# Merge control (ctrol_obj) and mutant (exp_obj) together (all):
all <- merge(ctrol_obj, y = exp_obj, add.cell.ids = c("CTRL", "MUT"), project = "scRNAseqI")
all@meta.data$sampleInfo <- all@meta.data$orig.ident
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^mt-")
head(all@meta.data)
table(all@meta.data$sampleInfo )

########################################
# check for basic info and data polish
########################################

#Make a Violin Plot for basic information: g1: total number of RNA molecules "nCount_RNA", g2: percent of mitochondrial DNA"percent.mt", 
# g3: total number of gene "nFeature_RNA":

VlnPlot(all, features = c( "nCount_RNA", "percent.mt","nFeature_RNA"), ncol = 3,group.by = "sampleInfo")
mito.drop <- isOutlier(all@meta.data$percent.mt, nmads=3, type="higher")
attributes(mito.drop)$thresholds
Gene.drop <- isOutlier(all@meta.data$nFeature_RNA, nmads=3, type="both")
attributes(Gene.drop)$thresholds
summary(all@meta.data$nFeature_RNA)
summary(all@meta.data$nFeature_RNA)
UMI.drop <- isOutlier(all@meta.data$nCount_RNA, nmads=3, type="both")
attributes(UMI.drop)$thresholds
UMI_lower=as.integer(attributes(UMI.drop)$thresholds[1])
UMI_upper=as.integer(attributes(UMI.drop)$thresholds[2])
Gene_lower=as.integer(attributes(Gene.drop)$thresholds[1])
Gene_upper=as.integer(attributes(Gene.drop)$thresholds[2])
mito_upper=round(attributes(mito.drop)$thresholds[2],3)

#Make a ggplot for basic information: g1: percent of mitochondrial DNA"percent.mt", g2: total number of gene "nFeature_RNA"
# g3: total number of RNA molecules "nCount_RNA":

g1= ggplot(all@meta.data, aes(x=percent.mt)) + geom_histogram(binwidth=.5)+
  geom_vline(xintercept = round(mito_upper[[1]],2), linetype="dashed")+
  geom_text(aes(x = round(mito_upper[[1]],2)+2, y = -10, label = as.character(round(mito_upper[[1]],2))))
g2= ggplot(all@meta.data, aes(x=nFeature_RNA)) + geom_histogram(binwidth=100)+
  geom_vline(xintercept = Gene_upper[[1]], linetype="dashed")+
  geom_text(aes(x = Gene_upper[[1]]+2, y = -10, label = as.character(Gene_upper[[1]])))
g3= ggplot(all@meta.data, aes(x=nCount_RNA)) + geom_histogram(binwidth=100)+
  geom_vline(xintercept = UMI_upper[[1]], linetype="dashed")+
  geom_text(aes(x = UMI_upper[[1]]+2, y = -10, label = as.character(UMI_upper[[1]])))
g4 = VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
g5=plot1 + plot2

# Save ggplots in pdf and name it 'basic info'
pdf('Basic_info.pdf')
print ( g1+g2+g3)
g4
g5
dev.off()

################################################
# RUN PCA DIMENSIONALITY REDUCTION -> 1st METHOD
################################################

# Clean the "all" file by removing cells with nFeature_RNA > 200 & nFeature_RNA < & percent.mt < 10:

all.clean <- subset(all, subset = nFeature_RNA > 200 & nFeature_RNA <  Gene_upper[[1]] & percent.mt < 10)
all.clean <- SCTransform(all.clean, vars.to.regress = "percent.mt", verbose = FALSE)
all.clean <- RunPCA(all.clean, verbose = FALSE)
all.clean <- RunUMAP(all.clean, dims = 1:30, verbose = FALSE)
all.clean <- FindNeighbors(all.clean, dims = 1:30, verbose = FALSE)

# Find clusters in "all" based on resolutions (0.1,0.2,0.3,0.4,0.6, 0.8, 1.2): 

all.clean <- FindClusters(
  object = all.clean,
  reduction.type = "pca",
  resolution = c(0.1,0.2,0.3,0.4,0.6, 0.8, 1.2),
  dims.use = 1:30,
  save.SNN = TRUE
)
head(all.clean[[]])

# Obtain the graph with all the clusters (g6):
g6=clustree(all.clean, prefix = "SCT_snn_res.")
g6

# Find selected clusters according to resolution:
res_val="res_0.2"
all.clean <- FindClusters(all.clean, verbose = FALSE,resolution = 0.2)
head(all.clean[[]])

# DimPlot (Dimensional Reduction Plot is a 2D scatter plot where each point is a cell and its position depend on the resolution)

g7=DimPlot(all.clean,label = F,split.by = 'sampleInfo')

# pdf 'PCA analysis'

pdf('PCA_analysis.pdf')
g6
g7
dev.off()

# RData dataset 

scrna=all.clean
save(scrna,file = paste('Seurat_merged.RData',sep=''))
rm(scrna)

####################################################
# RUN FASTMNN DIMENSIONALITY REDUCTION -> 2nd METHOD
####################################################

all.clean.2 <- subset(all, subset = nFeature_RNA > 200 & nFeature_RNA <  Gene_upper[[1]] & percent.mt < 10)
all.clean.2 <- NormalizeData(all.clean.2)
all.clean.2 <- FindVariableFeatures(all.clean.2)
all.clean.2 <- RunFastMNN(object.list = SplitObject(all.clean.2, split.by = "sampleInfo"))
all.clean.2 <- RunUMAP(all.clean.2, reduction = "mnn", dims = 1:30)
all.clean.2 <- FindNeighbors(all.clean.2, reduction = "mnn", dims = 1:30)

# Find all clusters according to resolution

all.clean.2 <- FindClusters(
  object = all.clean.2,
  reduction = "mnn",
  resolution = c(0.1,0.2,0.3,0.4,0.6, 0.8, 1.2),
  dims.use = 1:30,
  save.SNN = TRUE
)

head(all.clean.2[[]])
g8=clustree(all.clean.2, prefix = "RNA_snn_res.")

# find selected cluster (choose among clustree results according to resolution)

res_val="res.0.1"
all.clean.2 <- FindClusters(all.clean.2, verbose = FALSE,resolution = 0.1)

# DimPlot (g10 is split by CTRL and MUT)

# DimPlot for all cells (control + mutant)
g9=DimPlot(all.clean.2,label =T,cols = cbPalette)+labs(title = 'Seurat_clusters FastMnn [res.0.1]' )

#DimPlot split by group
g10=DimPlot(all.clean.2,label = F,split.by = 'sampleInfo',cols = cbPalette)+
  labs(title = 'Seurat_clusters FastMnn [res.0.1]' )

# pdf FASTMNN Analysis

pdf('FastMNN_analysis.pdf')
print(g8)
print(g9)
print(g10)
dev.off()


###########################################################################
# FIND MARKERS AMONG CLUSTERS (markers are differentially expressed genes)
###########################################################################

# Make a selection of genes you want to keep 
# min_ptc selectes genes that are expressed in at least 25% of cells in the cluster of interest
# logfc selects genes that are different in the cluster of interest compare to other clusters, the difference should be at least 0.5 (in log2 scale)

min_pct=0.25
logfc=0.5

Gene_markers.clusters <- FindAllMarkers(all.clean.2, only.pos = T, min.pct = min_pct, logfc.threshold = logfc)
table(Gene_markers.clusters$cluster)
head(Gene_markers.clusters)

# Make a cvs file with selected markers according to min_pct and logfc

write.csv(Gene_markers.clusters,paste('FastMNN_analysis.csv'))

all.clean.2 <- ScaleData(all.clean.2)
scrna=all.clean.2
save(scrna,file = paste('Seurat_fastamnn_','.processed.RData',sep=''))
rm(scrna)

###############################################
# HEATMAPS UNBIASED AND DOTPLOTs UN- AND BIASED
###############################################

load("~/Folder with data/Seurat_fastamnn_.processed.RData")
merge=scrna # 
rm(scrna)

# Consider only top 10 genes expressed in each cluster: 

markers_top10.clusters= Gene_markers.clusters%>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)


# pdf 'Heatmaps and DotPlots'

pdf(file = paste('Heatmaps_and_DotPlots.pdf',sep=''),height = 13,width = 20)
print(DotPlot(all.clean.2, features = targetGenes)+ RotatedAxis() + FontSize(x.tex = 8,y.tex = 10) )
print(DoHeatmap(merge, features = markers_top10.clusters$gene,size = 5) + 
        FontSize(x.tex = 0,y.tex = 7)+labs(title = 'seurat_clusters: top 5 Genes' ))
DotPlot(merge,features =unique(markers_top10.clusters$gene),split.by = 'sampleInfo' )+
  RotatedAxis()+labs(title = 'seurat_clusters top 10 Genes [res.0.1]' )+ FontSize(x.tex = 6,y.tex = 8)
DotPlot(merge, features = unique(markers_top10.clusters$gene), cols = c("blue", "red"),
        dot.scale = 8, split.by = 'sampleInfo') + 
  RotatedAxis() + FontSize(x.tex = 6,y.tex = 8)
dev.off()
DotPlot(merge, features = unique(markers_top10.clusters$gene), cols = c("blue", "red"),
        dot.scale = 8, split.by = 'sampleInfo') + 
  RotatedAxis() + FontSize(x.tex = 6,y.tex = 8)
