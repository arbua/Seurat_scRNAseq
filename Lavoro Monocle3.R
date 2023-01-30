# Author: Arboit Lorenzo
#Date: September, 29th 2021

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

setwd("~/First_data")
targetGenes=unique(c("Insert genes of interest as their HGNC code"))


load("~/Folder with data/Seurat_fastamnn_.processed.RData")
scrna

DefaultAssay(scrna) <- "RNA"
table(scrna@meta.data$sampleInfo)
ctrl=scrna[,scrna$sampleInfo=='CTRL']
mutant=scrna[,scrna$sampleInfo=='MUT']


##############
# ctrl
##############

ctrl.raw.data <-as(as.matrix(GetAssayData(ctrl, slot = "counts")),'sparseMatrix')
gene_names <- row.names(ctrl.raw.data)
gene_short_name <- gene_names
gene_annotation <- cbind(gene_short_name, gene_names)
rownames(gene_annotation) <- gene_names
gene_annotation <- as.data.frame(gene_annotation)
ctrl.cds <- new_cell_data_set(expression_data =  ctrl.raw.data,
                              cell_metadata = (ctrl@meta.data),
                              gene_metadata = gene_annotation)
ctrl.cds <- preprocess_cds(ctrl.cds, num_dim = 50)
ctrl.cds <- align_cds(ctrl.cds, alignment_group = "batch")

# Step 3: Reduce the dimensions using UMAP

ctrl.cds <- reduce_dimension(ctrl.cds)
coul <- colorRampPalette(brewer.pal(8, "RdYlBu"))(25)
c1= plot_cells(ctrl.cds,genes=c("Main genes of interest"),
           label_groups_by_cluster=FALSE, 
           label_leaves=FALSE,show_trajectory_graph = F,
           label_branch_points=FALSE)+
  scale_color_gradientn(colours = rainbow(5))
c2= plot_cells(ctrl.cds,
           genes=targetGenes,
           # color_cells_by = "cluster",
           label_roots   =F,
           label_groups_by_cluster=F,
           show_trajectory_graph=F,
           label_leaves=FALSE,
           label_branch_points=FALSE)+
  scale_color_gradientn(colours = rainbow(5))
ctrl.cds <- cluster_cells(ctrl.cds)

# plot_cells(ctrl.cds, color_cells_by = "partition")

ctrl.cds <- learn_graph(ctrl.cds)
ctrl.cds <- order_cells(ctrl.cds)
c3= plot_cells(ctrl.cds,
           # genes=targetGenes,
           color_cells_by = "pseudotime",
           label_roots   =F,
           label_groups_by_cluster=F,
           show_trajectory_graph=T,
           label_leaves=FALSE,
           label_branch_points=FALSE)
c4= plot_cells(ctrl.cds,
           # genes=targetGenes,
           color_cells_by = "cluster",
           label_groups_by_cluster=F,
           show_trajectory_graph=F,
           label_leaves=T,
           label_branch_points=FALSE)

scrna <- AddMetaData(
  object = scrna,
  metadata = ctrl.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "ctrl.pseudotime"
)

pdf('Monocle3_CTRL.pdf')
c1
c2
c3
c4
dev.off()

##############
# mutant
##############

mutant.raw.data <-as(as.matrix(GetAssayData(mutant, slot = "counts")),'sparseMatrix')
gene_names <- row.names(mutant.raw.data)
gene_short_name <- gene_names
gene_annotation <- cbind(gene_short_name, gene_names)
rownames(gene_annotation) <- gene_names
gene_annotation <- as.data.frame(gene_annotation)
mutant.cds <- new_cell_data_set(expression_data =  mutant.raw.data,
                                cell_metadata = (mutant@meta.data),
                                gene_metadata = gene_annotation)
mutant.cds <- preprocess_cds(mutant.cds, num_dim = 50)
mutant.cds <- align_cds(mutant.cds, alignment_group = "batch")

# Step 3: Reduce the dimensions using UMAP

mutant.cds <- reduce_dimension(mutant.cds)
m1= plot_cells(mutant.cds,genes=c("Main genes of interest"),
           label_groups_by_cluster=FALSE, 
           label_leaves=FALSE,show_trajectory_graph = F,
           label_branch_points=FALSE)+
  scale_color_gradientn(colours = rainbow(5))
m2= plot_cells(mutant.cds,
           genes=targetGenes,
           # color_cells_by = "cluster",
           label_roots   =F,
           label_groups_by_cluster=F,
           show_trajectory_graph=F,
           label_leaves=FALSE,
           label_branch_points=FALSE)+
  scale_color_gradientn(colours = rainbow(5))
mutant.cds <- cluster_cells(mutant.cds)

# plot_cells(mutant.cds, color_cells_by = "partition")

mutant.cds <- learn_graph(mutant.cds)
mutant.cds <- order_cells(mutant.cds)
m3= plot_cells(mutant.cds,
           # genes=targetGenes,
           color_cells_by = "pseudotime",
           label_roots   =F,
           label_groups_by_cluster=F,
           show_trajectory_graph=T,
           label_leaves=FALSE,
           label_branch_points=FALSE)
m4= plot_cells(mutant.cds,
           # genes=targetGenes,
           color_cells_by = "cluster",
           label_groups_by_cluster=F,
           show_trajectory_graph=F,
           label_leaves=T,
           label_branch_points=FALSE)

scrna <- AddMetaData(
  object = scrna,
  metadata = mutant.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "mutant.pseudotime"
)

pdf('Monocle3_MUT.pdf')
m1
m2
m3
m4
dev.off()

################
# FINAL PLOT PDF
################

l1= FeaturePlot(scrna, c("ctrl.pseudotime", "mutant.pseudotime"), pt.size = 0.1) & scale_color_viridis_c(option = "plasma")

l2= FeaturePlot(scrna, features=c("Main genes of interest"),
            pt.size = 0.1) #& scale_color_viridis_c()

pdf('Monocle3.pdf')
l1
l2
dev.off()
