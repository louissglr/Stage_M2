# Script: Single-cell RNA-seq Quality Control and Filtering with Seurat
#
# Description:
# This script performs quality control (QC) and filtering on single-cell RNA sequencing (scRNA-seq) data using Seurat.
# It loads Seurat objects from RDS files, applies thresholds to filter cells based on metrics such as total counts,
# number of features detected, mitochondrial content, gene-to-UMI ratio, and doublet scores.
# The script generates density plots and UMAP visualizations before and after filtering,
# runs dimensionality reduction, clustering, and saves filtered Seurat objects and QC plots.
#
# Input:
# - RDS files containing Seurat objects located in `input_dir`
#
# Output:
# - QC plots saved as PDF files in `output_pdf_dir`
# - Filtered and clustered Seurat objects saved as RDS files in `output_rds_dir`
#
library(Seurat)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(SingleCellExperiment)
library(scDblFinder)
library(biomaRt)
library(gprofiler2)

input_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/02_NormalizedData"
output_pdf_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/03_FilteredData/quality_controls_plot"
output_rds_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/03_FilteredData"

if (!dir.exists(output_pdf_dir)) dir.create(output_pdf_dir, recursive = TRUE)
if (!dir.exists(output_rds_dir)) dir.create(output_rds_dir, recursive = TRUE)

rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

plot_density_with_threshold <- function(data, thresholds = NULL, fill_color, title, label_positions = NULL) {
  data_frame <- data.frame(value = data)
  x_label <- deparse(substitute(data))
  x_label <- sub(".*\\$", "", x_label)
  
  plot <- ggplot(data = data_frame, aes(x = value)) +
    geom_density(fill = fill_color, alpha = 0.5) +
    ggtitle(title) +
    xlab(x_label) +
    ylab("Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (!is.null(thresholds)) {
    for (i in seq_along(thresholds)) {
      position <- if (!is.null(label_positions) && length(label_positions) >= i) label_positions[i] else 0
      plot <- plot +
        geom_vline(xintercept = thresholds[i], linetype = "dashed", color = "black", size = 1) +
        geom_text(aes(x = thresholds[i] + position, y = 0.0001, label = as.character(thresholds[i])), 
                  vjust = -1, size = 4, color = "black")
    }
  }
  
  return(plot)
}

plot_umap <- function(seurat_obj, group_by) {
  DimPlot(seurat_obj, reduction = "umap", group.by = group_by) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

generate_grid_plots <- function(plot1, plot2) {
  grid.arrange(plot1, plot2, ncol = 2)
}

for (file in rds_files) {
  sample_name <- gsub("^seurat_(.*)\\.rds$", "\\1", basename(file))
  seurat_obj <- readRDS(file)
  assign(sample_name, seurat_obj)
  
  seuil_count_max <- 7000
  seuil_count_min <- 50
  seuil_feature_min <- 200
  seuil_feature_max <- 5000
  seuil_mito_max <- 15
  seuil_log10_gene_umi <- 0.7
  seuil_doublet_score <- 0.6
  
  seurat_obj$seuil_count <- seurat_obj$nCount_RNA > seuil_count_min & seurat_obj$nCount_RNA < seuil_count_max
  seurat_obj$seuil_feature <- seurat_obj$nFeature_RNA > seuil_feature_min & seurat_obj$nFeature_RNA < seuil_feature_max
  seurat_obj$seuil_mito <- seurat_obj$percent.mt < seuil_mito_max
  seurat_obj$seuil_geneumi <- seurat_obj$log10GenesPerUMI > seuil_log10_gene_umi
  seurat_obj$seuil_doublet <- seurat_obj$doublets.score < seuil_doublet_score
  
  seurat_obj$seuil_final <- seurat_obj$seuil_count & seurat_obj$seuil_feature & seurat_obj$seuil_mito & seurat_obj$seuil_geneumi & seurat_obj$seuil_doublet
  percentage_fail <- round(mean(!seurat_obj$seuil_final) * 100, 2)
  
  resolution <- seurat_obj@commands$FindClusters$resolution 
  num_dims <- length(seurat_obj@commands$RunUMAP$dims)
  
  umap_plot_final <- plot_umap(seurat_obj, "seuil_final") +
    ggtitle(paste("Percentage of cells failing QC:", percentage_fail, "%"))
  
  seurat_filtered <- subset(seurat_obj, subset = seuil_final == TRUE)
  
  seurat_filtered <- RunPCA(seurat_filtered, verbose = FALSE)
  seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:30)

  seurat_filtered <- FindClusters(seurat_filtered, resolution = c(seq(0.1, 1, by = 0.1)), verbose = FALSE)

  seurat_filtered$seurat_clusters <- seurat_filtered$RNA_snn_res.0.4

  Idents(seurat_filtered) <- seurat_filtered$seurat_clusters
  seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:30)
  
  num_cells_after_qc <- ncol(seurat_filtered)
  
  umap_plot_clusters <- plot_umap(seurat_filtered, "seurat_clusters") +
    ggtitle(paste("Clustering after QC - res.0.4 - dim = 1:30 - n=", num_cells_after_qc, "\nDims:", num_dims, "- Resolution:", resolution))
  
  output_pdf <- file.path(output_pdf_dir, paste0("QC_", sample_name, ".pdf"))
  output_rds <- file.path(output_rds_dir, paste0("seurat_filtered_", sample_name, ".rds"))
  
  pdf(output_pdf, height = 8.5, width = 11)
  
  generate_grid_plots(
    plot_density_with_threshold(
      seurat_obj$nCount_RNA,
      thresholds = c(seuil_count_min, seuil_count_max),
      fill_color = "lightblue",
      title = "nCount_RNA density",
      label_positions = c(100, 100)
    ),
    plot_umap(seurat_obj, "seuil_count") +
      ggtitle(paste("nCount_RNA UMAP\nDims:", num_dims, "- Resolution:", resolution))
  )
  
  generate_grid_plots(
    plot_density_with_threshold(
      seurat_obj$nFeature_RNA,
      thresholds = c(seuil_feature_min, seuil_feature_max),
      fill_color = "pink",
      title = "nFeature_RNA density",
      label_positions = c(200, 200)
    ),
    plot_umap(seurat_obj, "seuil_feature") +
      ggtitle(paste("nFeature_RNA UMAP\nDims:", num_dims, "- Resolution:", resolution))
  )
  
  generate_grid_plots(
    plot_density_with_threshold(
      seurat_obj$percent.mt,
      thresholds = c(seuil_mito_max),
      fill_color = "lightgreen",
      title = "% Mito density",
      label_positions = c(0.1)
    ),
    plot_umap(seurat_obj, "seuil_mito") +
      ggtitle(paste("% Mito UMAP\nDims:", num_dims, "- Resolution:", resolution))
  )
  
  grid.arrange(
    plot_density_with_threshold(seurat_obj$percent.ribo, fill_color = "lightcoral", title = "% Ribosomique Density"),
    ncol = 1
  )

  generate_grid_plots(
    plot_density_with_threshold(
      seurat_obj$log10GenesPerUMI,
      thresholds = c(seuil_log10_gene_umi),
      fill_color = "lightgrey",
      title = "log10GenesPerUMI density",
      label_positions = c(0.02)
    ),
    plot_umap(seurat_obj, "seuil_geneumi") +
      ggtitle(paste("log10GenesPerUMI UMAP\nDims:", num_dims, "- Resolution:", resolution))
  )
  
  generate_grid_plots(
    plot_density_with_threshold(
      seurat_obj$doublets.score,
      thresholds = c(seuil_doublet_score),
      fill_color = "lightblue",
      title = "Doublet Score Density",
      label_positions = c(0.02)
    ),
    plot_umap(seurat_obj, "seuil_doublet") +
      ggtitle(paste("Doublet Score UMAP\nDims:", num_dims, "- Resolution:", resolution))
  )
  
  grid.arrange(umap_plot_final, ncol = 1)
  
  grid.arrange(umap_plot_clusters, ncol = 1)
  
  if ("TdTomato" %in% rownames(seurat_filtered)) {
    plot_tdTomato <- FeaturePlot(seurat_filtered, features = "TdTomato", reduction = "umap") +
      ggtitle(paste("UMAP Representation of TdTomato Expression\nDims:", num_dims, "- Resolution:", resolution))
    print(plot_tdTomato)
  } else {
    message("TdTomato non présent dans l'objet Seurat, figure non générée.")
  }
  
  if ("SV40" %in% rownames(seurat_filtered)) {
    plot_sv40 <- FeaturePlot(seurat_filtered, features = "SV40", reduction = "umap") +
      ggtitle(paste("UMAP Representation of SV40 Expression\nDims:", num_dims, "- Resolution:", resolution))
    print(plot_sv40)
  } else {
    message("SV40 non présent dans l'objet Seurat, figure non générée.")
  }

  dev.off()
  saveRDS(seurat_filtered, file = output_rds)
}
