# Script: Generate Loupe Browser Files from Seurat Objects using loupeR
#
# Description:
# This script processes multiple merged single-cell RNA-seq Seurat objects to create Loupe Browser-compatible files.
# It reads Seurat objects from a directory, assigns unique barcodes by combining sample identifiers with original barcodes,
# then uses the `loupeR` package to generate `.cloupe` files for visualization in the 10x Genomics Loupe Browser.
#
# Input:
# - Reference Seurat object with barcode information (`ref_file`)
# - Multiple Seurat RDS files located in `input_dir`
#
# Output:
# - Loupe Browser files saved to `output_dir`
#
library(loupeR)

input_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData"
ref_file <- "C:/Users/louis/Desktop/Stage/scrnaseq/output/loupeR/ex_10x.rds"
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/05_CreateloupeR"

seurat_obj <- readRDS(ref_file)
barcodes_obj <- colnames(seurat_obj)
n_obj <- length(barcodes_obj)

rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

for (rds_path in rds_files) {
  seurat <- readRDS(rds_path)
  DefaultAssay(seurat) <- "RNA"
  seurat$original_barcode <- colnames(seurat)
  seurat$file_source <- basename(rds_path)
  
  n_cells <- ncol(seurat)
  recycled_barcodes <- rep(barcodes_obj, length.out = n_cells)
  new_barcodes <- paste0(seurat$orig.ident, "_", recycled_barcodes)
  colnames(seurat) <- new_barcodes
  
  file_base <- tools::file_path_sans_ext(basename(rds_path))
  output_name <- paste0("cogaps_", file_base)
  
  setwd(output_dir)
  create_loupe_from_seurat(seurat, output_name = output_name)
}

