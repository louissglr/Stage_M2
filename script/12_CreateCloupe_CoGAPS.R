# Script: Export Seurat Objects with CoGAPS Results to Loupe Browser Format
#
# Description:
# This script converts Seurat objects containing CoGAPS pattern information into `.cloupe`-compatible
# files using the `loupeR` package. These outputs can be loaded into the 10x Genomics Loupe Browser
# for interactive exploration of CoGAPS-based clustering and dimensionality reduction.
#
# Inputs:
# - A directory of Seurat RDS files with CoGAPS results (`AddCoGAPSToSeurat` output).
# - A reference Seurat object (`ex_10x.rds`) to extract and recycle barcodes.
#
# Outputs:
# - One Loupe-compatible dataset per input Seurat object, written to the `output_dir` as `.cloupe` files.
#
# Key Steps:
# - Read reference barcodes from a template Seurat object.
# - Iterate over CoGAPS-augmented Seurat objects.
# - Assign recycled barcodes and generate unique cell IDs.
# - Export each Seurat object as a Loupe-compatible file using `create_loupe_from_seurat()`.
#
library(loupeR)
library(Seurat)

input_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat"
ref_file <- "C:/Users/louis/Desktop/Stage/scrnaseq/output/loupeR/ex_10x.rds"
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/12_CreateCloupe_CoGAPS"

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

