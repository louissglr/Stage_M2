# Script: Barplot and Cell Count Summary by CoGAPS Patterns and Time Point
#
# Description:
# This script analyzes the distribution of cells across CoGAPS-derived patterns and experimental 
# time points (`orig.ident`) within a Seurat object. It visualizes the results using a stacked barplot 
# and exports the summarized data for further interpretation.
#
# Objectives:
# - Visualize the number of cells associated with each CoGAPS pattern across time points.
# - Label patterns only when they represent ≥ 10% of the cells at a given time point.
# - Export the raw table of cell counts per pattern and time point.
#
# Inputs:
# - `seurat_cogaps_res_10p.rds`: A Seurat object that includes CoGAPS pattern assignments (`pattern_cogaps`)
#   and time point metadata (`orig.ident`).
#
# Outputs:
# - `nb_cellules_par_patterns_et_orig.ident.csv`: A CSV file listing the number of cells per pattern and time point.
#
library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)

seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")
seurat_obj$pattern_cogaps <- gsub("Pattern_", "", seurat_obj$pattern_cogaps)

df <- seurat_obj@meta.data %>%
  count(pattern_cogaps, orig.ident) %>%
  group_by(orig.ident) %>%
  mutate(
    prop = n / sum(n),
    label = ifelse(prop >= 0.1, as.character(pattern_cogaps), NA)  
  ) %>%
  arrange(orig.ident, desc(pattern_cogaps)) %>%
  mutate(pos = cumsum(n) - n / 2) %>%
  ungroup()

barplot <- ggplot(df, aes(x = factor(orig.ident), y = n, fill = factor(pattern_cogaps))) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(y = pos, label = label), color = "white", size = 5, fontface = "bold", na.rm = TRUE) + 
  labs(
    title = "Number of Cells by Patterns and Time Point (orig.ident)",
    x = "Time Point (orig.ident)",
    y = "Number of Cells",
    fill = "Patterns"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(color = "white", size = 16, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "white", size = 16), 
    axis.title.x = element_text(color = "white", size = 18, face = "bold"),  
    axis.title.y = element_text(color = "white", size = 18, face = "bold"),  
    legend.text = element_text(color = "white", size = 14),
    legend.title = element_text(color = "white", size = 14, face = "bold"),
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"),
    plot.title = element_text(hjust = 0.5, color = "white", size = 20, face = "bold") 
  )

print(barplot)

ggsave(
  "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation/barplot_cogaps.png", 
  plot = barplot, 
  width = 10,    
  height = 8,   
  dpi = 300     
)


df <- seurat_obj@meta.data %>%
  count(orig.ident, pattern_cogaps) %>%
  rename(
    patterns = pattern_cogaps,
    nb_cellules = n
  ) %>%
  arrange(orig.ident, patterns)

print(head(df))
write_csv(df, "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation/nb_cellules_par_patterns_et_orig.ident.csv")
