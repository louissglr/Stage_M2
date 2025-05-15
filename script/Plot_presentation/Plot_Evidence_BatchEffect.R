# Script: Barplot Visualization of Cluster and Pattern Composition in Seurat Object
#
# Description:
# This script generates barplots to visualize the composition of:
#   1. Seurat clusters across CoGAPS patterns 2 and 6.
#   2. Time points (orig.ident) within Seurat clusters 2, 3, and 9.
#   3. Time points within clusters 12, 10, and 4 (data assumed to be in `df_cluster_12_10_4`).
#
# Inputs:
# - `seurat_cogaps_res_10p.rds`: A Seurat object with CoGAPS results and metadata.
# - Metadata columns used: `pattern_cogaps`, `seurat_clusters`, and `orig.ident`.
#
# Outputs:
# - Three horizontal barplots saved as PNG files in:
#     * `barplot_cogaps.png`
#     * `barplot_cluster_orig.png`
#     * `barplot_cluster_12_10_4.png` (uses `df_cluster_12_10_4`, assumed precomputed)
#
library(ggplot2)
library(dplyr)
library(scales)

cogaps <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_10p.rds")

cogaps$pattern_cogaps <- gsub("Pattern_", "", cogaps$pattern_cogaps)

df <- cogaps@meta.data %>%
  filter(pattern_cogaps %in% c("2", "6")) %>%    
  count(pattern_cogaps, seurat_clusters) %>%
  group_by(pattern_cogaps) %>%
  mutate(
    prop = n / sum(n),
    label = ifelse(prop >= 0.05, as.character(seurat_clusters), NA)  #label si prop >= 5%
  ) %>%
  ungroup()

df <- df %>%
  group_by(pattern_cogaps) %>%
  arrange(pattern_cogaps, desc(seurat_clusters)) %>%
  mutate(pos = cumsum(prop) - prop / 2) %>%
  ungroup()

barplot <- ggplot(df, aes(x = factor(pattern_cogaps), y = prop, fill = seurat_clusters)) +
  geom_bar(stat = "identity", color = "white") +
  geom_text(aes(y = pos, label = label), color = "white", size = 6, fontface = "bold", na.rm = TRUE) +  # Texte dans barres plus grand
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Proportion of Seurat Clusters for Patterns 2 and 6",
    x = "Dominant Pattern",
    y = "Proportion of Seurat Clusters"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(color = "white", size = 16),              
    axis.text.y = element_text(color = "white", size = 18, face = "bold"),
    axis.title.x = element_text(color = "white", size = 20),          
    axis.title.y = element_text(color = "white", size = 20),            
    legend.text = element_text(color = "white", size = 14),
    legend.title = element_text(color = "white", size = 14),
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"),
    plot.title = element_text(hjust = 0.5, color = "white", size = 16, face = "bold")
  )

print(barplot)

ggsave(
  filename = "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation/barplot_cogaps.png",
  plot = barplot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "black"
)

df_cluster_orig <- cogaps@meta.data %>%
  filter(seurat_clusters %in% c("2", "3", "9")) %>%
  count(seurat_clusters, orig.ident) %>%
  group_by(seurat_clusters) %>%
  mutate(
    prop = n / sum(n),
    label = ifelse(prop >= 0.1, as.character(orig.ident), NA)  
  ) %>%
  ungroup()

df_cluster_orig <- df_cluster_orig %>%
  group_by(seurat_clusters) %>%
  arrange(seurat_clusters, desc(orig.ident)) %>%
  mutate(pos = cumsum(prop) - prop / 2) %>%
  ungroup()

barplot_cluster_orig <- ggplot(df_cluster_orig, aes(x = factor(seurat_clusters), y = prop, fill = orig.ident)) +
  geom_bar(stat = "identity", color = "white") +
  geom_text(aes(y = pos, label = label), color = "white", size = 6, fontface = "bold", na.rm = TRUE) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Proportion of orig.ident in Seurat Clusters 2, 3, and 9",
    x = "Seurat Cluster",
    y = "Proportion of orig.ident",
    fill = "Time Point (orig.ident)"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(color = "white", size = 16),
    axis.text.y = element_text(color = "white", size = 18, face = "bold"),
    axis.title.x = element_text(color = "white", size = 20),
    axis.title.y = element_text(color = "white", size = 20),
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
    plot.title = element_text(hjust = 0.5, color = "white", size = 18, face = "bold")
  )

print(barplot_cluster_orig)

ggsave(
  filename = "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation/barplot_cluster_orig.png",
  plot = barplot_cluster_orig,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "black"
)

barplot_cluster_12_10_4 <- ggplot(df_cluster_12_10_4, aes(x = factor(seurat_clusters), y = prop, fill = orig.ident)) +
  geom_bar(stat = "identity", color = "white") +
  geom_text(aes(y = pos, label = label), color = "white", size = 6, fontface = "bold", na.rm = TRUE) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Proportion of orig.ident in Clusters 12, 10 and 4",
    x = "Seurat Cluster",
    y = "Proportion of orig.ident",
    fill = "Time Point (orig.ident)"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(color = "white", size = 16),
    axis.text.y = element_text(color = "white", size = 18, face = "bold"),
    axis.title.x = element_text(color = "white", size = 20),
    axis.title.y = element_text(color = "white", size = 20),
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
    plot.title = element_text(hjust = 0.5, color = "white", size = 24, face = "bold")  # Titre centré et plus gros
  )
