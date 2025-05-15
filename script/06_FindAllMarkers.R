# Script: Identification and Export of Cluster Marker Genes at Multiple Resolutions with Seurat
#
# Description:
# This script reads merged Seurat objects and performs marker gene identification for each clustering resolution.
# For each resolution, it finds all positive marker genes for clusters using FindAllMarkers,
# separates significant markers (adjusted p-value < 0.05), and exports results to Excel files.
# Each resolution has two sheets per file: all markers and significant markers.
#
# Input:
# - Merged Seurat RDS files located in the specified directory
#
# Output:
# - Excel files containing marker gene lists for each clustering resolution, saved in the output directory
#
library(Seurat)
library(openxlsx)
library(dplyr)

files <- list.files("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData", pattern = "*.rds", full.names = TRUE)

resolutions <- seq(0.1, 1, by = 0.1)

for (file in files) {
  
  obj <- readRDS(file)
  obj <- JoinLayers(obj)
  
  wb <- createWorkbook()
  
  for (res in resolutions) {
    res_label <- paste0("RNA_snn_res.", res)
    
    if (res_label %in% colnames(obj@meta.data)) {
      Idents(obj) <- obj[[res_label]][, 1]
      
      #markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      markers <- FindAllMarkers(obj, only.pos = TRUE)
      markers <- markers %>%
        dplyr::select(gene, everything()) 
      
      markers_significant <- markers %>% filter(p_val_adj < 0.05)
      
      addWorksheet(wb, sheetName = paste0("res_", res, "_significant"))
      writeData(wb, sheet = paste0("res_", res, "_significant"), x = markers_significant)
      
      markers_all <- markers  
      
      addWorksheet(wb, sheetName = paste0("res_", res, "_all"))
      writeData(wb, sheet = paste0("res_", res, "_all"), x = markers_all)
    }
  }
  
  output_file <- paste0("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/06_FindAllMarkers/Markers_", basename(file), ".xlsx")
  
  saveWorkbook(wb, file = output_file, overwrite = TRUE)
  
  print(paste("Fichier Excel créé pour", basename(file)))
}
