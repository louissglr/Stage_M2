# Script: Single-cell RNA-seq preprocessing with Seurat and scDblFinder
# 
# Description:
# This script processes single-cell RNA sequencing (scRNA-seq) data using Seurat.
# It loads Seurat objects, detects doublets, calculates mitochondrial and ribosomal content,
# performs normalization (NormalizeData) and clustering, and assigns cell cycle scores.
#
# Input:
# - RDS files containing Seurat objects stored in `input_dir`
#
# Output:
# - Processed Seurat objects with doublet classification, mitochondrial and ribosomal percentages,
#   cell cycle scoring, and clustering, saved as RDS files in `output_rds_dir`

# Load required libraries
library(Seurat)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(SingleCellExperiment)
library(scDblFinder)
library(biomaRt)
library(gprofiler2)

# Define input and output directories
input_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/01_CreateObject"
output_rds_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/02_NormalizedData"

# Create output directory if it does not exist
if (!dir.exists(output_rds_dir)) dir.create(output_rds_dir, recursive = TRUE)

# List all RDS files in the input directory
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

# Process each RDS file
for (file in rds_files) {
  sample_name <- gsub("^seurat_(.*)\\.rds$", "\\1", basename(file))
  seurat_obj <- readRDS(file)
  assign(sample_name, seurat_obj)
  
  # Convert Seurat object to SingleCellExperiment
  sceobj <- as.SingleCellExperiment(seurat_obj)
  
  # Predict doublets using scDblFinder
  doublet_res <- scDblFinder(sceobj)
  seurat_obj$doublets.class <- doublet_res$scDblFinder.class
  seurat_obj$doublets.score <- doublet_res$scDblFinder.score
  
  # Compute mitochondrial gene percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  # Compute ribosomal gene percentage
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")  

  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  # Normalize the data using LogNormalization
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identify highly variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # Scale data while regressing out unwanted variation
  seurat_obj <- ScaleData(seurat_obj)
  
  # Perform cell cycle scoring
  cc_seurat <- Seurat::cc.genes.updated.2019
  mmus_s = gorth(cc_seurat$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
  mmus_g2m = gorth(cc_seurat$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
  
  seurat_obj <- CellCycleScoring(
    object = seurat_obj,
    s.features = mmus_s,
    g2m.features = mmus_g2m,
    assay = 'RNA',
    seed = 123  
  )
  
  # Equivalent to Min-cells when creating Seurat object
  counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  genes_to_keep <- rownames(counts_matrix)[Matrix::rowSums(counts_matrix > 0) >= 3]
  seurat_obj <- subset(seurat_obj, features = genes_to_keep)
  
  # Re-scale the data, now including cell cycle scores
  #seurat_obj <- ScaleData(seurat_obj)
  
  # Perform PCA and clustering
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  
  # Save the processed Seurat object
  output_rds <- file.path(output_rds_dir, paste0("seurat_", sample_name, ".rds"))
  saveRDS(seurat_obj, file = output_rds)
}
