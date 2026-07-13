# data from human 14.000 PBMC w. 2000 genes/features from two conditions: before+After interferon treatment
# interferon-treatment/stimulated PBMC and control PBMCcells(not interferon-stimulated) from the same person

####################################
# 1. Imports ----
####################################
options(timeout = 10000)
library(Seurat)
library(SeuratData)
library(patchwork)
library(tidyverse)
library(sctransform)

# restart session if errors
# SeuratData::InstallData("ifnb.SeuratData", timeout = 10000)

# load dataset
ifnb <- LoadData("ifnb")

UpdateSeuratObject(ifnb)


####################################
# 2. Preprocessing----
####################################

ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
# run standard preprocessing workflow
ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- RunUMAP(ifnb, dims = 1:30)

####################################
# 3. Intermediate Processing w/o Integration-----
####################################

ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")

# create a dim. reduction w. name unintegrated_clusters
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")

####################################
# 4. Unintegrated Analysis ----
####################################

# problem:
# interferons have strong effect on gene-expression of nearly all pbmc celltypes
# clustering w/o integration is useless
DimPlot(ifnb, reduction = "umap", group.by = c("seurat_annotations", "unintegrated_clusters"))

DimPlot(ifnb, reduction = "umap", group.by = c("seurat_annotations", "unintegrated_clusters", "stim"))


####################################
# 5. Integrated Analysis ----
####################################

ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method = "SCT")

# process integrated data
ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")

# result: treatment now has no influence on the UMPAP clustering 
DimPlot(ifnb, reduction = "umap", group.by = c("seurat_annotations", "seurat_clusters", "stim"))





# Prep step to FindMarkers w. SC-Tr / perform differential expression
ifnb <- PrepSCTFindMarkers(ifnb) # Function to prepare FindMarkers

# assign list to new metadata col celltype.stim 
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_") # create the interaction identities

# assign the interaction identities
Idents(ifnb) <- "celltype.stim" 

DimPlot(ifnb, reduction = "umap", split.by = "stim")

# FindMarker genes that are expressed as cd14 celltype's response to interferons, dh:
# FindMarkers Genes that separate "CD14 Mono_STIM" from "CD14 Mono_CTRL" d.h.
CD14Mono_stim_markers <- FindMarkers(ifnb, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL")

# show differential expression(diffex) of top 15 marker for cd14 stim vs ctrl
# high diffex: log2fc > 1 and <-1
# small diffex log2fc  between -1 and 1 or look in table of findconservedmarkers
head(CD14Mono_stim_markers, n = 15)

# save interferone response markergenes to a vector/list
CD14Mono_markers_top15 <- rownames(head(CD14Mono_stim_markers, n = 15))

head(CD14Mono_markers_top15, 15)

# create pseudobulks
aggregate_ifnb <- AggregateExpression(ifnb, group.by = c("seurat_annotations", "stim"), return.seurat = TRUE)

# p1_u an p3_u are helperplot thats adds the dots/genes to plt p2_u
p1_u <- CellScatter(aggregate_ifnb, "CD14 Mono_CTRL", "CD14 Mono_STIM", highlight = CD14Mono_markers_top15)
# this adds the labels
p2_u <- LabelPoints(plot = p1_u, points = CD14Mono_markers_top15, repel = TRUE, max.overlaps = 15)

p3_u <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", highlight = CD14Mono_markers_top15)
p4_u <- LabelPoints(plot = p3_u, points = CD14Mono_markers_top15, repel = TRUE, max.overlaps = 15)

# highlight the CD14 markers of the CD14 interferon response
# left: compare CD14 before and after interferon stimulation
# right: see how theese markers are affected in the CD4 interferon response
p2_u + p4_u

# unbiased gene labelling - cd14 markers on two monocyte cell types
p1_u <- CellScatter(aggregate_ifnb, "CD14 Mono_CTRL", "CD14 Mono_STIM", highlight = CD14Mono_markers_top15)
p2_u <- LabelPoints(plot = p1_u, points = CD14Mono_markers_top15, repel = TRUE, max.overlaps = 15)
p3_u <- CellScatter(aggregate_ifnb, "CD16 Mono_CTRL", "CD16 Mono_STIM", highlight = CD14Mono_markers_top15)
p4_u <- LabelPoints(plot = p3_u, points = CD14Mono_markers_top15, repel = TRUE, max.overlaps = 15)
p2_u+p4_u

# get the top 15 marker sets for visualisation
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM", ident.2 = "B_CTRL")
head(b.interferon.response, n = 15)
b_top15 <- rownames(head(b.interferon.response, n = 15))

# unbiased gene labelling - B markers on both B-cell types
p1_u <- CellScatter(aggregate_ifnb, "B_CTRL", "B_STIM", highlight = b_top15)
p2_u <- LabelPoints(plot = p1_u, points = b_top15, repel = TRUE, max.overlaps = 15)
p3_u <- CellScatter(aggregate_ifnb, "B Activated_CTRL", "B Activated_STIM", highlight = b_top15)
p4_u <- LabelPoints(plot = p3_u, points = b_top15, repel = TRUE, max.overlaps = 15)
p2_u+p4_u

# Get the marker sets for visualisation
cd4_naive.interferon.response <- FindMarkers(ifnb, ident.1 = "CD4 Naive T_STIM", ident.2 = "CD4 Naive T_CTRL")
head(cd4_naive.interferon.response, n = 15)
cd4_naive_top15 <- rownames(head(cd4_naive.interferon.response, n = 15))
cd4_memory.interferon.response <- FindMarkers(ifnb, ident.1 = "CD4 Memory T_STIM", ident.2 = "CD4 Memory T_CTRL")
head(cd4_memory.interferon.response, n = 15)
cd4_memory_top15 <- rownames(head(cd4_memory.interferon.response, n = 15))

# unbiased gene labelling - cd4 naive vs cd4 memory interferon responses
p1_u <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", highlight = cd4_naive_top15)
p2_u <- LabelPoints(plot = p1_u, points = cd4_naive_top15, repel = TRUE, max.overlaps = 15)
p3_u <- CellScatter(aggregate_ifnb, "CD4 Memory T_CTRL", "CD4 Memory T_STIM", highlight = cd4_naive_top15)
p4_u <- LabelPoints(plot = p3_u, points = cd4_naive_top15, repel = TRUE, max.overlaps = 15)
p2_u+p4_u

# unbiased gene labelling - cd4 naive vs cd8 interferon responses
p1_u <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", highlight = cd4_naive_top15)
p2_u <- LabelPoints(plot = p1_u, points = cd4_naive_top15, repel = TRUE, max.overlaps = 15)
p3_u <- CellScatter(aggregate_ifnb, "CD8 T_CTRL", "CD8 T_STIM", highlight = cd4_naive_top15)
p4_u <- LabelPoints(plot = p3_u, points = cd4_naive_top15, repel = TRUE, max.overlaps = 15)
p2_u+p4_u

# unbiased gene labelling - cd4 memory vs cd8 interferon responses
p1_u <- CellScatter(aggregate_ifnb, "CD4 Memory T_CTRL", "CD4 Memory T_STIM", highlight = cd4_memory_top15)
p2_u <- LabelPoints(plot = p1_u, points = cd4_memory_top15, repel = TRUE, max.overlaps = 15)
p3_u <- CellScatter(aggregate_ifnb, "CD8 T_CTRL", "CD8 T_STIM", highlight = cd4_memory_top15)
p4_u <- LabelPoints(plot = p3_u, points = cd4_memory_top15, repel = TRUE, max.overlaps = 15)
p2_u+p4_u

# Genes chosen to show pattern recognition:
# CD3D and GNLY are canonical cell type markers (for T cells and NK/CD8 T cells) 
# But They are not really affected by interferon stim
# IFI6 and ISG15 are affected by interferon / response genes in all cell types
FeaturePlot(ifnb, features = c("CD3D", "GNLY", "IFI6"), 
            split.by = "stim", max.cutoff = 3, 
            cols = c("grey",  "red"), 
            reduction = "umap")
          
# CD14 (the CD14/16 cell marker gene), is downreg. 
# bacteria invades, produces toxinens, cd14 is a detector for toxines, once the toxines are detected the inflammation starts, this downregulates cd14 / to balance the process/ to prevent over Inflammation
# CXCL10 (CD16 marker, upreg) is upreguted in several CELL-TYPES 
FeaturePlot(ifnb, features = c("CD14", "ISG15", "CXCL10"), 
            split.by = "stim", max.cutoff = 3, 
            cols = c("grey",  "red"), 
            reduction = "umap")
          
DimPlot(ifnb, reduction = "umap", group.by="seurat_annotations")

# exactly same genes and the same conclusion as violin plots:
#    CD3D and GNLY are canonical cell type markers (for T cells and NK/CD8 T cells) 
#    But They are not really affected by interferon stim
plots <- VlnPlot(ifnb, features = c("CD3D", "GNLY"), split.by = "stim", group.by = "seurat_annotations",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# IFI6 and ISG15 are affected by interferon / response genes in all cell types
plots <- VlnPlot(ifnb, features = c("IFI6", "ISG15"), split.by = "stim", group.by = "seurat_annotations",
                 pt.size = 0, combine = FALSE) 
wrap_plots(plots = plots, ncol = 1)

# CD14 (the CD14/16 cell marker gene), is downreg 
# bacteria invades, produces toxinens, cd14 is a detector for toxines, once the toxines are detected the inflammation starts, this downregulates cd14 / to balance the process/ to prevent over Inflammation
# CXCL10 (CD16 marker, upreg) is upreguted in several CELL-TYPES 
plots <- VlnPlot(ifnb, features = c("CD14", "CXCL10"), split.by = "stim", group.by = "seurat_annotations",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

sink(paste0("session_info_Integration", Sys.Date(), ".txt"))
sessionInfo()
sink()

