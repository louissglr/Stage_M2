# Script: Visualization of ImSig Module Scores by Seurat Clusters
#
# Description:
# This script loads a combined Seurat object and a mouse immune signature gene list (ImSig).
# It calculates a module score for the ImSig signature in each cell,
# then creates a violin plot of the ImSig scores grouped by Seurat clusters.
# The plot uses a custom black background theme and includes cell count annotations per cluster.
# The final plot is saved as a PNG image.
#
# Input:
# - Combined Seurat object RDS file
# - ImSig mouse signature gene list (text file)
#
# Output:
# - PNG file with violin plot of ImSig scores by cluster
#
library(dplyr)
library(ggplot2)
library(patchwork)

seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds")
signature_imsig <- read.table("C:/Users/louis/Desktop/Stage/scrnaseq/signatures/ImSig2_mouse.txt", header = FALSE, stringsAsFactors = FALSE)$V1
seurat_obj <- AddModuleScore(seurat_obj, features = list(signature_imsig), name = "ImSig_Score")

theme_black_background <- theme_minimal() +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"),
    axis.text = element_text(color = "white", size = 14),
    axis.title = element_text(color = "white"),
    panel.grid = element_blank(),
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    plot.title = element_text(color = "white", hjust = 0.5)
  )

add_cell_count <- function(p, data, group_by) {
  counts <- data %>%
    group_by(.data[[group_by]]) %>%
    summarise(
      count = n(),
      max_val = max(ImSig_Score1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_pos = max_val + 0.05)  
  
  p + geom_text(
    data = counts,
    aes(x = .data[[group_by]], y = y_pos, label = paste("n =", count)),
    color = "white",
    size = 2.7,
    fontface = "bold",
    inherit.aes = FALSE
  )
}

p1 <- VlnPlot(
  object = seurat_obj,
  features = "ImSig_Score1",
  group.by = "seurat_clusters",
  pt.size = 0
) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
  labs(title = "ImSig Score by Seurat Cluster") +
  theme_black_background +
  theme(legend.position = "none")

max_y1 <- max(seurat_obj@meta.data$ImSig_Score1[!is.na(seurat_obj@meta.data$seurat_clusters)], na.rm = TRUE)
p1 <- p1 + scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(NA, max_y1 * 1.15))

p1 <- add_cell_count(p1, seurat_obj@meta.data, "seurat_clusters")

ggsave("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/07_SignatureImSigCluster/ImSig_ViolinPlot_SeuratClusters.png", plot = p1, width = 6, height = 6, dpi = 300, bg = "black")
