# Charger les deux objets
seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")
a <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_2p.rds")

# Assurer que les cellules sont dans le même ordre
stopifnot(all(colnames(seurat_obj) == colnames(a)))

# Copier la colonne pattern depuis a vers seurat_obj
seurat_obj$pattern_cogaps <- a$pattern_cogaps

# Nettoyage du nom (enlève "Pattern_")
seurat_obj$pattern_cogaps2 <- gsub("Pattern_", "", seurat_obj$pattern_cogaps)

# Définir les identités pour Seurat
Idents(seurat_obj) <- "pattern_cogaps2"

# Assay par défaut si nécessaire
DefaultAssay(seurat_obj) <- "CoGAPS"

# Visualisation
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
seurat_obj$pattern_cogaps <- gsub("Pattern_", "", seurat_obj$pattern_cogaps)
p <- DimPlot(seurat_obj, reduction = "umap", group.by = "pattern_cogaps", label = TRUE, label.size = 8, label.color = "white") +
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

print(p)
umap_output_dir <- "C:/Users/louis/Desktop"
ggsave(
  filename = file.path(umap_output_dir, "cogaps_umap_patterns.png"),
  plot = p,
  bg = "black",
  width = 7,
  height = 6,
  dpi = 300
)
