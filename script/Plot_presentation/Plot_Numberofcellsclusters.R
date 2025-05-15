# Script: Barplot and Cell Count Summary by Cluster and Time Point
#
# Description:
# This script analyzes and visualizes the distribution of cells across Seurat clusters and
# time points (`orig.ident`) using a stacked barplot. It also exports a summary table 
# containing the number of cells per cluster and time point.
#
# Objectives:
# - Visualize the number of cells from each cluster across different experimental conditions (orig.ident).
# - Add labels to clusters representing at least 10% of a given time point’s composition.
# - Export a table of raw cell counts for further downstream analysis or reporting.
#
# Inputs:
# - `seurat_combined.rds`: A merged Seurat object containing `seurat_clusters` and `orig.ident` metadata.
#
# Outputs:
# - `barplot.png`: A high-contrast stacked barplot showing cell counts per time point and cluster.
# - `nb_cellules_par_cluster_et_orig.ident.csv`: A CSV file summarizing cell counts by cluster and time point.
#
library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)
seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds")


df <- seurat_obj@meta.data %>%
  count(seurat_clusters, orig.ident) %>%
  group_by(orig.ident) %>%
  mutate(
    prop = n / sum(n),
    label = ifelse(prop >= 0.1, as.character(seurat_clusters), NA)  
  ) %>%
  arrange(orig.ident, desc(seurat_clusters)) %>%
  mutate(pos = cumsum(n) - n / 2) %>%
  ungroup()


barplot <- ggplot(df, aes(x = factor(orig.ident), y = n, fill = factor(seurat_clusters))) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(y = pos, label = label), color = "white", size = 5, fontface = "bold", na.rm = TRUE) +  # Increased label size
  labs(
    title = "Number of Cells by Cluster and Time Point (orig.ident)",
    x = "Time Point (orig.ident)",
    y = "Number of Cells",
    fill = "Cluster"
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
  "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation/barplot.png", 
  plot = barplot, 
  width = 10,   
  height = 8,  
  dpi = 300      
)

df <- seurat_obj@meta.data %>%
  count(orig.ident, seurat_clusters) %>%
  rename(
    clusters = seurat_clusters,
    nb_cellules = n
  ) %>%
  arrange(orig.ident, clusters)

print(head(df))
write_csv(df, "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation/nb_cellules_par_cluster_et_orig.ident.csv")
