library(Seurat)
library(ggplot2)

# Lire l'objet Seurat combiné
seurat <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")
seurat$pattern_cogaps <- gsub("Pattern_", "", seurat$pattern_cogaps)
# Créer le répertoire de sortie si nécessaire
umap_output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation"
if (!dir.exists(umap_output_dir)) {
  dir.create(umap_output_dir, recursive = TRUE)
}

# Générer le UMAP avec fond noir et axes blancs, par samples
p <- DimPlot(seurat, reduction = "umap", group.by = "orig.ident") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    
    panel.grid = element_blank(),
    
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

# Sauvegarder le graphique
ggsave(
  filename = file.path(umap_output_dir, "cogaps_umap_origident.png"),
  plot = p,
  bg = "black",
  width = 7,
  height = 6,
  dpi = 300
)

p <- DimPlot(seurat, reduction = "umap", group.by = "pattern_cogaps", label = TRUE, label.size = 8, label.color = "white") +
  
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    
    panel.grid = element_blank(),
    
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )


# Sauvegarder le graphique
ggsave(
  filename = file.path(umap_output_dir, "cogaps_umap_patterns.png"),
  plot = p,
  bg = "black",
  width = 7,
  height = 6,
  dpi = 300
)

