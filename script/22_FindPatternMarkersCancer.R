# Script: Extract CoGAPS Pattern Markers and Save to Excel Workbook
#
# Description:
# This script iterates over multiple CoGAPS result directories to extract pattern-specific gene markers
# using the `patternMarkers` function. Each set of markers corresponding to a pattern is organized
# into an Excel sheet, with one sheet per CoGAPS result directory.
#
# Inputs:
# - A root directory (`res_dir`) containing subdirectories with CoGAPS results (`cogaps.rds` files).
#
# Outputs:
# - A single Excel file (`pattern_markers.xlsx`) containing one worksheet per subdirectory.
#   Each sheet includes the list of genes identified as markers for each pattern.
#
# Key Steps:
# - Sort subdirectories using natural (mixed) ordering.
# - For each subdirectory:
#   - Load the `cogaps.rds` result.
#   - Extract gene markers per pattern using `patternMarkers`.
#   - Normalize data frame length across patterns for formatting.
#   - Write the results to a worksheet in a shared Excel file.
#
############### Import libraries ################################################
library(CoGAPS)
library(openxlsx)
library(gtools)  # Pour mixedsort()

############### Path setup ######################################################
res_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/20_CoGAPS_cancer/"
output_file <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/22_FindPatternMarkersCancer/pattern_markers.xlsx"

wb <- createWorkbook()

subdirs <- list.dirs(res_dir, recursive = FALSE)

subdirs_sorted <- mixedsort(subdirs)

for (subdir in subdirs_sorted) {
  cogaps_file <- file.path(subdir, "cogaps.rds")
  
  if (file.exists(cogaps_file)) {
    message("Traitement de : ", cogaps_file)
    
    cogaps <- readRDS(cogaps_file)
    
    patternMarkerResults <- patternMarkers(cogaps, threshold = "all")
    
    df <- data.frame(lapply(patternMarkerResults$PatternMarkers, function(x) {
      length(x) <- max(lengths(patternMarkerResults$PatternMarkers))
      return(x)
    }), stringsAsFactors = FALSE)
    
    colnames(df) <- paste0("Pattern_", seq_along(df))
    
    sheet_name <- substr(basename(subdir), 1, 31)
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, x = df)
    
  } else {
    warning("Fichier non trouvé : ", cogaps_file)
  }
}

saveWorkbook(wb, output_file, overwrite = TRUE)
message("Fichier Excel sauvegardé : ", output_file)
