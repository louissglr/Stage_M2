# Script: Integrating CoGAPS Patterns into Seurat Object and UMAP Visualization
#
# Description:
# This script loads a previously merged Seurat object and integrates CoGAPS results
# (learned pattern matrices) into the object. For each CoGAPS result stored in a subdirectory,
# it adds the patterns as a new assay ("CoGAPS"), sets it as the default, computes the dominant
# pattern per cell, adds it to metadata, and constructs a dimensional reduction object.
# A UMAP is then computed and visualized for each CoGAPS run, grouped by dominant pattern.
#
# Inputs:
# - A merged Seurat object containing single-cell RNA-seq data.
# - A directory containing multiple subdirectories, each with a `cogaps.rds` result.
#
# Outputs:
# - Updated Seurat objects with integrated CoGAPS results saved as RDS files.
# - UMAP plots saved as PDF files, visualizing CoGAPS pattern clustering.
#
###############Import libraries###################################################
library(Seurat)
library(CoGAPS)
library(dplyr)
library(ggplot2)
library(forcats)
library(knitr)
library(ComplexHeatmap)
library(grid)
library(openxlsx)

###############Path setup#########################################################
# Set paths
res_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/08_CoGAPS/"
seurat_file <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds"
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/"

# Load Seurat object
seurat <- readRDS(seurat_file)

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

###############Loop through subdirectories and process CoGAPS#####################
# Get list of subdirectories
subdirs <- list.dirs(res_dir, recursive = FALSE)

for (subdir in subdirs) {
  cogaps_file <- file.path(subdir, "cogaps.rds")
  
  # Check if cogaps.rds exists in the subdirectory
  if (file.exists(cogaps_file)) {
    message("Processing: ", cogaps_file)
    
    # Load CoGAPS result
    cogaps <- readRDS(cogaps_file)
    
    # Extract patterns and make sure order matches
    patterns_in_order <- t(cogaps@sampleFactors[colnames(seurat), ])
    
    # Add CoGAPS assay to a copy of Seurat object
    seurat_temp <- seurat
    seurat_temp[["CoGAPS"]] <- CreateAssayObject(counts = patterns_in_order)
    DefaultAssay(seurat_temp) <- "CoGAPS"
    
    # Identify dominant pattern per cell
    dominant_patterns <- apply(patterns_in_order, 2, function(x) {
      rownames(patterns_in_order)[which.max(x)]
    })
    
    # Add to metadata
    seurat_temp@meta.data$pattern_cogaps <- dominant_patterns[colnames(seurat_temp)]
    seurat_temp[["cogaps"]] <- CreateDimReducObject(
      embeddings = t(as.matrix(patterns_in_order)), 
      key = "GAPS_",
      assay = "CoGAPS"
    )
    
    # Run UMAP on CoGAPS reduction
    seurat_temp <- RunUMAP(seurat_temp, reduction = "cogaps", dims = 1:nrow(patterns_in_order))
    
    folder_name <- basename(subdir)
    output_file <- file.path(output_dir, paste0("seurat_cogaps_", folder_name, ".rds"))
    saveRDS(seurat_temp, output_file)
    message("Saved to: ", output_file)
    
    pdf_file <- file.path(output_dir, paste0("umap_cogaps_", folder_name, ".pdf"))
    pdf(pdf_file, width = 7, height = 6)
    print(
      DimPlot(seurat_temp, reduction = "umap", group.by = "pattern_cogaps") +
        ggtitle(paste("UMAP - CoGAPS patterns:", folder_name)) +
        theme_minimal()
    )
    dev.off()
    message("UMAP PDF saved to: ", pdf_file)
    
  } else {
    warning("No cogaps.rds found in: ", subdir)
  }
}
