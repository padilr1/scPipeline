#!/usr/bin/env Rscript
#' @title scAnalyze
#'
#' @description Seurat workflow for preprocessing, quality control, clustering, and cell-type annotation of scRNA-seq data
#' @param mat a character path to a tab-delimited counts matrix (genes in column 1; cells with corresponding read counts in remaining columns).
#' @param proj_name a character string specifying the project name
#' @param seed an integer seed value for reproducibility. Default 42.
#' @param min_features a numeric to keep cells with nFeature_RNA (number of detected genes) > `min_features`. Default 100. This is a typical minimum for Drop-seq; 10x Chromium datasets may use a different cutoff
#' @param max_features a numeric to keep cells with nFeature_RNA (number of detected genes) < `max_features`. Default 6000. This is a typical maximum for Drop-seq; 10x Chromium datasets may use a different cutoff
#' @param max_counts a numeric to keep cells with nCount_RNA (number of UMIs) < `max_counts`. Default 10000. This is a typical maximum for Drop-seq; 10x Chromium datasets may use a different cutoff
#' @param max_percent_mito a numeric between 0 and 1.0 to keep cells with percent.mito < `max_percent_mito`. Default 0.10.
#' @param max_percent_crystal a numeric between 0 and 1.0 to keep cells with percent.crystal < `max_percent_crystal`. Default 0.025.
#' @param pcs_use an integer vector specifying which principal components to use. Default 1:20
#' @param res a numeric controlling clustering granularity. Higher values result in more fine-grained clusters; lower values yield broader groups. Default 0.5
#' @return a list with elements: Seurat object, markers and plots.
#'
#' @examples
#' \dontrun{
#' scAnalyze(mat = "counts_matrix.txt", proj_name = "example_scRNAseq_analysis", pcs_use = 1:6, res = 0.5)
#' }
#' @details developed and tested using Seurat v5.3.0.

scAnalyze <- function(mat,
                      proj_name,
                      seed = 42,
                      min_features = 100,
                      max_features = 6000,
                      max_counts = 10000,
                      max_percent_mito = 0.10,
                      max_percent_crystal = 0.025,
                      pcs_use = 1:20,
                      res = 0.5) {
  # ---- Dependencies needed
  needed <- c("data.table", "Matrix", "Seurat", "patchwork", "dplyr", "tibble", "ggplot2")
  miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "))

  # ---- Seurat version check
  if (utils::packageVersion("Seurat") < "5.3.0") {
    warning("Tested with Seurat >= 5.3.0; your version is ", utils::packageVersion("Seurat"))
  }

  # ---- Input checks
  if (!is.character(mat) || length(mat) != 1L || !file.exists(mat)) {
    stop("`mat` must be a path to an existing tab-delimited file.")
  }
  if (!is.character(proj_name) || length(proj_name) != 1L) {
    stop("`proj_name` must be a single character string.")
  }

  # Enforce numeric types
  if (!is.numeric(min_features)) stop("`min_features` must be numeric.")
  if (!is.numeric(max_features)) stop("`max_features` must be numeric.")
  if (!is.numeric(max_counts)) stop("`max_counts` must be numeric.")
  if (!is.numeric(max_percent_mito)) stop("`max_percent_mito` must be numeric.")
  if (!is.numeric(max_percent_crystal)) stop("`max_percent_crystal` must be numeric.")
  if (!is.numeric(res)) stop("`res` must be numeric.")

  if (!is.numeric(pcs_use) || any(pcs_use %% 1 != 0) || any(pcs_use <= 0)) {
    stop("`pcs_use` must be a vector of positive integers (e.g., 1:20).")
  }

  # ensure values are positive
  if (min_features < 0) stop("`min_features` must be >= 0.")
  if (max_features <= 0) stop("`max_features` must be > 0.")
  if (max_counts <= 0) stop("`max_counts` must be > 0.")
  if (res <= 0) stop("`res` must be > 0.")

  if (max_percent_mito < 0 || max_percent_mito > 1) {
    stop("`max_percent_mito` must be in [0, 1].")
  }
  if (max_percent_crystal < 0 || max_percent_crystal > 1) {
    stop("`max_percent_crystal` must be in [0, 1].")
  }
  if (min_features >= max_features) {
    stop("`min_features` must be less than `max_features`.")
  }

  # --- set seed for reproducibility
  set.seed(seed)

  # --- load counts (genes in column 1)
  DGEmatrix <- data.table::fread(mat, sep = "\t") %>%
    data.frame(row.names = 1)
  # stop if number of columns is less than 2
  if (ncol(DGEmatrix) < 2) {
    stop("Input matrix must have ≥2 columns: gene identifiers and read counts.")
  }

  # --- create Seurat object, discarding cells that have less than 100 genes
  # and discarding genes that are expressed in less than 3 cells.
  Seurat_object <- Seurat::CreateSeuratObject(DGEmatrix,
    assay = "RNA",
    min.cells = 3,
    min.features = 100,
    project = proj_name
  )

  # --- QC metrics (mito/crystallin)
  # pull the counts matrix
  counts <- Seurat::GetAssayData(Seurat_object, assay = "RNA", slot = "counts")

  # get a list of mitochondrial genes: all genes starting with 'MT'
  mito.genes <- grep("^MT-", rownames(counts), value = TRUE)

  # compute library size
  libsize <- Matrix::colSums(counts)

  # compute proportions - measurement of mitochondrial gene proportion
  percent.mito <- if (length(mito.genes)) {
    Matrix::colSums(counts[mito.genes, , drop = FALSE]) / libsize
  } else {
    rep(0, ncol(counts))
  }

  # add the information back to the seurat object
  Seurat_object <- Seurat::AddMetaData(
    object = Seurat_object, metadata = percent.mito,
    col.name = "percent.mito"
  )

  # get a list of crystal genes (all genes starting by 'CRY') to remove crystalin contamination
  crystal.genes <- grep("^CRY[AB]", rownames(x = Seurat_object@assays$RNA), value = TRUE)

  # compute proportions
  percent.crystal <- if (length(crystal.genes)) {
    Matrix::colSums(counts[crystal.genes, , drop = FALSE]) / libsize
  } else {
    rep(0, ncol(counts))
  }

  # add the information back to the seurat object
  Seurat_object <- Seurat::AddMetaData(
    object = Seurat_object, metadata = percent.crystal,
    col.name = "percent.crystal"
  )

  # --- QC plots prior to filtering (pre)
  # # of genes, # of UMIs, and mitochondrial/crystal proportion
  fts <- c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal")
  preFilt_vlnPlot <- Seurat::VlnPlot(object = Seurat_object, ncol = 4, pt.size = 0, features = fts)

  # --- Summary stats for the sample prior to filtering (pre):
  summary_stats_before_filtering <- tibble::tibble(
    total_cells  = nrow(Seurat_object@meta.data),
    mean_n_genes = mean(Seurat_object@meta.data$nFeature_RNA),
    sd_n_genes   = sd(Seurat_object@meta.data$nFeature_RNA),
    max_n_genes  = max(Seurat_object@meta.data$nFeature_RNA),
    min_n_genes  = min(Seurat_object@meta.data$nFeature_RNA),
    mean_UMI     = mean(Seurat_object@meta.data$nCount_RNA),
    sd_UMI       = sd(Seurat_object@meta.data$nCount_RNA),
    max_UMI      = max(Seurat_object@meta.data$nCount_RNA),
    min_UMI      = min(Seurat_object@meta.data$nCount_RNA)
  ) %>% dplyr::mutate_all(function(x) round(x, 2))

  # --- Filtering
  Seurat_object <- subset(Seurat_object, subset = nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    nCount_RNA < max_counts &
    percent.mito < max_percent_mito &
    percent.crystal < max_percent_crystal)

  postFilt_vlnPlot <- Seurat::VlnPlot(object = Seurat_object, pt.size = 0, ncol = 4, features = fts)

  # --- Summary stats for the sample post-filtering:
  summary_stats_after_filtering <- tibble::tibble(
    total_cells  = nrow(Seurat_object@meta.data),
    mean_n_genes = mean(Seurat_object@meta.data$nFeature_RNA),
    sd_n_genes   = sd(Seurat_object@meta.data$nFeature_RNA),
    max_n_genes  = max(Seurat_object@meta.data$nFeature_RNA),
    min_n_genes  = min(Seurat_object@meta.data$nFeature_RNA),
    mean_UMI     = mean(Seurat_object@meta.data$nCount_RNA),
    sd_UMI       = sd(Seurat_object@meta.data$nCount_RNA),
    max_UMI      = max(Seurat_object@meta.data$nCount_RNA),
    min_UMI      = min(Seurat_object@meta.data$nCount_RNA)
  ) %>% dplyr::mutate_all(function(x) round(x, 2))

  print(
    paste0(
      "Number of filtered cells: ",
      summary_stats_before_filtering$total_cells -
        summary_stats_after_filtering$total_cells
    )
  )

  # --- Normalize dataset
  Seurat_object <- Seurat::NormalizeData(Seurat_object)

  # before proceeding, we note markers for cell cycle
  Seurat_object <- Seurat::CellCycleScoring(Seurat_object,
    set.ident = FALSE,
    s.features = cc.genes$s.genes,
    g2m.features = cc.genes$g2m.genes
  )

  # --- SCTransform: transforms the data and regresses out unwanted variation
  vs <- c("nFeature_RNA", "percent.mito", "S.Score", "G2M.Score")
  Seurat_object <- Seurat::SCTransform(Seurat_object,
    verbose = T,
    vars.to.regress = vs
  )

  # --- Post-SCTransform QC scatter
  p1 <- Seurat::FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mito")
  p2 <- Seurat::FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  postSC_plt <- (p1 | p2) + ggplot2::theme(legend.position = "none")

  # prepare seurat object for VariableFeaturePlot
  Seurat_object[["SCT"]]@meta.features <- Seurat::SCTResults(Seurat_object[["SCT"]], slot = "feature.attributes")[, c("gmean", "variance", "residual_variance")]
  Seurat_object[["SCT"]]@meta.features$variable <- F
  Seurat_object[["SCT"]]@meta.features[VariableFeatures(Seurat_object[["SCT"]]), "variable"] <- F
  colnames(Seurat_object[["SCT"]]@meta.features) <- paste0("sct.", colnames(Seurat_object[["SCT"]]@meta.features))

  # --- Label and plot the top 10 most variable features
  top10_variableFeatures <- Seurat::VariableFeaturePlot(Seurat_object, selection.method = "sct", assay = "SCT") %>%
    LabelPoints(points = head(VariableFeatures(Seurat_object), 10), repel = T) &
    theme(legend.position = "none")

  # --- PCA
  Seurat_object <- Seurat::RunPCA(Seurat_object, pcs.compute = 100, do.print = F)

  # genes correlated with PCs 1 & 2
  correlatedGenes_PC1and2 <- Seurat::VizDimLoadings(Seurat_object, dims = 1:2, reduction = "pca")

  # alternative representation of genes highly correlated with PC1
  heatmap_correlatedGenes_PC1 <- Seurat::DimHeatmap(Seurat_object, dims = pcs_use, cells = 500, balanced = T)

  # --- Elbowplot to inspect the standard deviation of each PC
  elbowPlot_PC <- Seurat::ElbowPlot(Seurat_object)

  # assess effect of cell cycle
  PCA_cellCycle <- Seurat::DimPlot(Seurat_object, reduction = "pca", group.by = "Phase")

  # --- Non-linear approaches: UMAP and tSNE
  Seurat_object <- Seurat::RunUMAP(Seurat_object, dims = pcs_use)
  Seurat_object <- Seurat::RunTSNE(Seurat_object, dims = pcs_use)

  # cluster the cells based on UMAP reduction
  Seurat_object <- Seurat::FindNeighbors(Seurat_object,
    dims = 1:2,
    reduction = "umap"
  )
  # the resolution can be adjusted to tweak clustering accuracy
  Seurat_object <- Seurat::FindClusters(Seurat_object,
    resolution = res,
    reduction = "umap"
  )

  # plot UMAP
  dimPlotUMAP <- Seurat::DimPlot(Seurat_object, reduction = "umap", label = T)

  # plot tSNE
  dimPlotTSNE <- Seurat::DimPlot(Seurat_object, reduction = "tsne", label = T)

  # plot metrics
  vlnPlotMetrics <- Seurat::VlnPlot(object = Seurat_object, pt.size = 0.1, ncol = 1, features = fts)


  # --- DE markers (switch to RNA for interpretability)
  # apply differential expression analysis to find marker genes higher expressed in every given cluster as compared to all remaining cells
  # then infer the cell type based on current knowledge on those genes
  Seurat::DefaultAssay(Seurat_object) <- "RNA"
  Seurat_object_markers <- Seurat::FindAllMarkers(Seurat_object,
    only.pos = TRUE,
    min.pct = 0.25, logfc.threshold = 0.25
  )

  # get the top 10 markers in each cluster
  top10_markers <- Seurat_object_markers %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 10) %>%
    ungroup()

  dotPlot_top10markers <- Seurat::DotPlot(
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

  # --- Compile plots
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

  # --- Workflow outputs
  out <- list(
    seurat = Seurat_object,
    markers = Seurat_object_markers,
    plots = list(QC = QC_plots, linear = linear_plots, nonlinear = nonlinear_plots)
  )

  return(out)
}
