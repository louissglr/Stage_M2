# Script: Annotate CoGAPS Patterns with ImSig Signature Scores and Visualize on UMAP
#
# Description:
# This script evaluates immune-related activity across CoGAPS-derived patterns using the ImSig mouse
# gene signature. It adds module scores for the ImSig signature to a CoGAPS-annotated Seurat object
# and visualizes the results using UMAP plots.
#
# Inputs:
# - `seurat_cogaps_res_10p.rds`: A Seurat object with CoGAPS patterns and UMAP embedding.
# - `ImSig2_mouse.txt`: A text file containing the ImSig mouse gene signature (one line, comma-separated).
#
# Outputs:
# - UMAP plots showing:
#   - Cells grouped by dominant CoGAPS pattern (`pattern_cogaps`).
#   - Cells grouped by sample (`orig.ident`).
#   - ImSig module score projected on CoGAPS UMAP embedding.
#
# Key Steps:
# - Load the Seurat object with CoGAPS results.
# - Score cells based on the ImSig gene set using `AddModuleScore()`.
# - Set the CoGAPS assay and clustering identity.
# - Visualize UMAPs with `DimPlot()` and `FeaturePlot()`.
###############Import library###################################################
library(Seurat)
library(CoGAPS)
library(dplyr)
library(ggplot2)
library(forcats)
library(knitr)
library(ComplexHeatmap)
library(grid)
library(openxlsx)
cogaps <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_10p.rds")
signature_imsig <- read.table("C:/Users/louis/Desktop/Stage/scrnaseq/signatures/ImSig2_mouse.txt", header = FALSE, stringsAsFactors = FALSE)$V1
DefaultAssay(cogaps) <- "RNA" 
cogaps <- AddModuleScore(
  cogaps, 
  features = list(signature_imsig), 
  name = "ImSig_Score"
)
DefaultAssay(cogaps) <- "CoGAPS"
Idents(cogaps) <- "pattern_cogaps"
DimPlot(cogaps, reduction = "umap", group.by = "pattern_cogaps", label=T)
DimPlot(cogaps, reduction = "umap", group.by = "orig.ident", label=T)
FeaturePlot(
  cogaps,
  features = "ImSig_Score1",
  reduction = "umap",  # UMAP calculé à partir de CoGAPS
  coord.fixed = TRUE,
  label=T
)
