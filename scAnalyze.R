#!/usr/bin/env Rscript
#' Title
#' @title Pre-processing, quality control, identifying cell populations and assigning identities of single-cell RNA samples using Seurat
#'
#' @description Using Seurat for processing, quality control, identifying cell populations and differential expression analysis.
#' @param proj_name a character string specifying the project name
#' @param mat a matrix of read counts per gene, per cell
#'
#' @return Seurat object and analysis plots
#' @export
#'
#'
#' @examples
#' \dontrun{
#' scAnalyze(mat = "", proj_name = "")
#' }
library(data.table)
library(tidyverse)
library(Seurat)
library(patchwork)
library(reactable)

scAnalyze <- function(mat, proj_name) {
  set.seed(42)
  DGEmatrix <- data.table::fread(mat, sep = "\t") %>%
    data.frame(row.names = 1)

  # create our base Seurat object, discarding cells that have less than 100 genes
  # and discarding genes that are expressed in less than 3 cells.
  Seurat_object <- Seurat::CreateSeuratObject(DGEmatrix,
    assay = "RNA",
    min.cells = 3, min.features = 100,
    project = proj_name
  )

  # pull the counts matrix
  counts <- Seurat::GetAssayData(Seurat_object, assay = "RNA", slot = "counts")

  # get a list of mitochondrial genes: all genes starting with 'MT'
  mito.genes <- grep("^MT-", rownames(counts), value = TRUE)

  # compute proportions - measurement of mitochondrial gene proportion
  percent.mito <- Matrix::colSums(counts[mito.genes, ]) /
    Matrix::colSums(counts)

  # add the information back to the seurat object
  Seurat_object <- Seurat::AddMetaData(
    object = Seurat_object, metadata = percent.mito,
    col.name = "percent.mito"
  )

  # get a list of crystal genes (all genes starting by 'CRY') to remove crystalin contamination
  crystal.genes <- grep("^CRY[AB]", rownames(x = Seurat_object@assays$RNA), value = TRUE)

  # compute proportions
  percent.crystal <- Matrix::colSums(counts[crystal.genes, ]) /
    Matrix::colSums(counts)

  # add the information back to the seurat object
  Seurat_object <- Seurat::AddMetaData(
    object = Seurat_object, metadata = percent.crystal,
    col.name = "percent.crystal"
  )

  # display distribution of some metrics:
  # # of genes, # of UMIs, and mitochondrial/crystal proportion
  fts <- c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal")
  preFilt_vlnPlot <- Seurat::VlnPlot(object = Seurat_object, ncol = 4, pt.size = 0, features = fts)

  # Get some summary stats for the sample prior to filtering:
  summary_stats_before_filtering <- tibble(
    total_cells  = nrow(Seurat_object@meta.data),
    mean_n_genes = mean(Seurat_object@meta.data$nFeature_RNA),
    sd_n_genes   = sd(Seurat_object@meta.data$nFeature_RNA),
    max_n_genes  = max(Seurat_object@meta.data$nFeature_RNA),
    min_n_genes  = min(Seurat_object@meta.data$nFeature_RNA),
    mean_UMI     = mean(Seurat_object@meta.data$nCount_RNA),
    sd_UMI       = sd(Seurat_object@meta.data$nCount_RNA),
    max_UMI      = max(Seurat_object@meta.data$nCount_RNA),
    min_UMI      = min(Seurat_object@meta.data$nCount_RNA)
  ) %>% mutate_all(function(x) round(x, 2))

  # We'll define some cutoffs to remove cells that deviate too much from the rest of the
  # cells in the sample. These can be hard cut-offs or dataset-specific ones based on SD.
  # Sometimes those represent multiplets, or just abnormal cells like crystalin cells
  Seurat_object <- subset(Seurat_object, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 &
    nCount_RNA < 10000 & percent.mito < 0.10 & percent.crystal < 0.025)
  postFilt_vlnPlot <- Seurat::VlnPlot(object = Seurat_object, pt.size = 0, ncol = 4, features = fts)

  # Get some summary stats for the sample post-filtering:
  summary_stats_after_filtering <- tibble(
    total_cells  = nrow(Seurat_object@meta.data),
    mean_n_genes = mean(Seurat_object@meta.data$nFeature_RNA),
    sd_n_genes   = sd(Seurat_object@meta.data$nFeature_RNA),
    max_n_genes  = max(Seurat_object@meta.data$nFeature_RNA),
    min_n_genes  = min(Seurat_object@meta.data$nFeature_RNA),
    mean_UMI     = mean(Seurat_object@meta.data$nCount_RNA),
    sd_UMI       = sd(Seurat_object@meta.data$nCount_RNA),
    max_UMI      = max(Seurat_object@meta.data$nCount_RNA),
    min_UMI      = min(Seurat_object@meta.data$nCount_RNA)
  ) %>% mutate_all(function(x) round(x, 2))

  print(
    paste0(
      "Number of filtered cells: ",
      summary_stats_before_filtering$total_cells -
        summary_stats_after_filtering$total_cells
    )
  )

  # normalize dataset
  Seurat_object <- Seurat::NormalizeData(Seurat_object)

  # before proceeding, we note markers for cell cycle
  Seurat_object <- Seurat::CellCycleScoring(Seurat_object,
    set.ident = FALSE,
    s.features = cc.genes$s.genes,
    g2m.features = cc.genes$g2m.genes
  )

  # we then go on to apply a transforms the data and regresses out unwanted variation
  vs <- c("nFeature_RNA", "percent.mito", "S.Score", "G2M.Score")
  Seurat_object <- Seurat::SCTransform(Seurat_object,
    verbose = T,
    vars.to.regress = vs
  )

  postSC_plt <- Seurat::FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mito") +
    FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") &
    theme(legend.position = "none")

  # a bit of wrangling to prepare for VariableFeaturePlot
  Seurat_object[["SCT"]]@meta.features <- SCTResults(Seurat_object[["SCT"]], slot = "feature.attributes")[, c("gmean", "variance", "residual_variance")]
  Seurat_object[["SCT"]]@meta.features$variable <- F
  Seurat_object[["SCT"]]@meta.features[VariableFeatures(Seurat_object[["SCT"]]), "variable"] <- F
  colnames(Seurat_object[["SCT"]]@meta.features) <- paste0("sct.", colnames(Seurat_object[["SCT"]]@meta.features))

  # label top 10 most variable features
  top10_variableFeatures <- Seurat::VariableFeaturePlot(Seurat_object, selection.method = "sct", assay = "SCT") %>%
    LabelPoints(points = head(VariableFeatures(Seurat_object), 10), repel = T) &
    theme(legend.position = "none")

  # perform PCA
  Seurat_object <- RunPCA(Seurat_object, pcs.compute = 100, do.print = F)

  # genes correlated with PCs 1 & 2
  correlatedGenes_PC1and2 <- VizDimLoadings(Seurat_object, dims = 1:2, reduction = "pca")

  # alternative representation of genes highly correlated with PC1
  heatmap_correlatedGenes_PC1 <- DimHeatmap(Seurat_object, dims = 1:15, cells = 500, balanced = T)

  # inspect the standard deviation of each PC
  elbowPlot_PC <- ElbowPlot(Seurat_object)

  # assess effect of cell cycle
  PCA_cellCycle <- DimPlot(Seurat_object, reduction = "pca", group.by = "Phase")

  # non-linear approaches: UMAP and tSNE
  Seurat_object <- RunUMAP(Seurat_object, dims = 1:20)
  Seurat_object <- RunTSNE(Seurat_object, dims = 1:20)

  # cluster the cells based on UMAP reduction
  Seurat_object <- FindNeighbors(Seurat_object,
    dims = 1:2,
    reduction = "umap"
  )
  # the resolution can be adjusted to tweak clustering accuracy
  Seurat_object <- FindClusters(Seurat_object,
    resolution = 0.5,
    reduction = "umap"
  )

  # plot UMAP
  dimPlotUMAP <- DimPlot(Seurat_object, reduction = "umap", label = T)

  # plot tSNE
  dimPlotTSNE <- DimPlot(Seurat_object, reduction = "tsne", label = T)

  # plot metrics
  vlnPlotMetrics <- VlnPlot(object = Seurat_object, pt.size = 0.1, ncol = 1, features = fts)

  # identifying cell types
  # apply differential expression analysis to find marker genes higher expressed in every given cluster as compared to all remaining cells
  # then infer the cell type based on current knowledge on those genes.
  Seurat_object.markers <- FindAllMarkers(Seurat_object,
    only.pos = TRUE,
    min.pct = 0.25, logfc.threshold = 0.25
  )

  # get the top 10 markers in each cluster
  top10_markers <- Seurat_object.markers %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 10) %>%
    ungroup()

  dotPlot_top10markers <- DotPlot(
    Seurat_object,
    assay = NULL,
    unique(top10_markers$gene),
    cols = c("blue", "red"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    group.by = NULL,
    split.by = NULL,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

  # compile plots
  QC_plots <- list(
    preFilt_vlnPlot = preFilt_vlnPlot,
    postFilt_vlnPlot = postFilt_vlnPlot,
    postSC_plt = postSC_plt
  )

  linear_plots <- list(
    top10_variableFeatures = top10_variableFeatures,
    correlatedGenes_PC1and2 = correlatedGenes_PC1and2,
    heatmap_correlatedGenes_PC1 = heatmap_correlatedGenes_PC1,
    elbowPlot_PC = elbowPlot_PC,
    PCA_cellCycle = PCA_cellCycle
  )

  nonlinear_plots <- list(
    dimPlotUMAP = dimPlotUMAP,
    dimPlotTSNE = dimPlotTSNE,
    vlnPlotMetrics = vlnPlotMetrics,
    dotPlot_top10markers = dotPlot_top10markers
  )

  # output
  out <- list(
    seurat = Seurat_object,
    plots = list(QC = QC_plots, linear = linear_plots, nonlinear = nonlinear_plots)
  )

  return(out)
}
