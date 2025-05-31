# Script: DotPlot Visualization of Selected Genes in Seurat Object
#
# Description:
# This script generates customized DotPlots to visualize expression of selected genes 
# (e.g., SV40, TdTomato, Krt14, and Vim) across cell clusters from a merged Seurat object.
# The DotPlots use a dark theme for enhanced contrast and publication-style visuals.
#
# Inputs:
# - `seurat_combined.rds`: A merged Seurat object containing RNA assay data and clustering information.
#
# Genes Visualized:
# - First plot: SV40, TdTomato, Krt14
# - Second plot: SV40, TdTomato, Krt14, Vim
#
# Outputs:
# - Two DotPlots displayed in the viewer, colored by average expression and sized by percent expression.
library(Seurat)
library(ggplot2)

seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds")

output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

##### Premier plot : SV40, TdTomato, Krt14
features1 <- c("SV40", "TdTomato", "Krt14")

p1 <- DotPlot(seurat_obj, features = features1) + RotatedAxis()

p1 <- p1 +
  scale_color_gradient(low = "blue", high = "red") + 
  labs(
    x = "Genes",
    y = "Clusters",
    color = "Average expression",
    size = "Percent expressed",
    title = "DotPlot - SV40, TdTomato, Krt14"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "white", size = 16, face = "bold"),
    axis.text.y = element_text(color = "white", size = 18, face = "bold"),
    axis.title.x = element_text(color = "white", size = 18),
    axis.title.y = element_text(color = "white", size = 18),
    plot.title = element_text(hjust = 0.5, color = "white", size = 16, face = "bold"),
    legend.text = element_text(color = "white", size = 10),
    legend.title = element_text(color = "white", size = 11),
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white")
  )

p1$layers[[1]]$aes_params$stroke <- 0.5  
p1$layers[[1]]$aes_params$color <- "white"  

# Sauvegarde plot 1
ggsave(filename = file.path(output_dir, "DotPlot_SV40_TdTomato_Krt14.png"), plot = p1,
       width = 10, height = 8, dpi = 300)

##### Deuxième plot : SV40, TdTomato, Krt14, Vim
features2 <- c("SV40", "TdTomato", "Krt14", "Vim")

p2 <- DotPlot(seurat_obj, features = features2) + RotatedAxis()

p2 <- p2 +
  scale_color_gradient(low = "blue", high = "red") +  
  labs(
    x = "Genes",
    y = "Clusters",
    color = "Average expression",
    size = "Percent expressed",
    title = "DotPlot - SV40, TdTomato, Krt14, Vim"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "white", size = 16, face = "bold"),
    axis.text.y = element_text(color = "white", size = 18, face = "bold"),
    axis.title.x = element_text(color = "white", size = 18),
    axis.title.y = element_text(color = "white", size = 18),
    plot.title = element_text(hjust = 0.5, color = "white", size = 16, face = "bold"),
    legend.text = element_text(color = "white", size = 10),
    legend.title = element_text(color = "white", size = 11),
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white")
  )

p2$layers[[1]]$aes_params$stroke <- 0.5  
p2$layers[[1]]$aes_params$color <- "white" 

# Sauvegarde plot 2
ggsave(filename = file.path(output_dir, "DotPlot_SV40_TdTomato_Krt14_Vim.png"), plot = p2,
       width = 10, height = 8, dpi = 300)
