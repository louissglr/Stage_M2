library(Seurat)
library(ggplot2)

seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_10p.rds")
DefaultAssay(seurat_obj) <- "RNA" 
seurat_obj$pattern_cogaps <- gsub("Pattern_", "", seurat_obj$pattern_cogaps)
Idents(seurat_obj) <- "pattern_cogaps"

# Gènes à visualiser
features <- c("SV40", "TdTomato", "Krt14")

# Générer le DotPlot standard
p <- DotPlot(seurat_obj, features = features) + RotatedAxis()

# Modifier uniquement le thème et les contours sans changer les échelles
p <- p +
  scale_color_gradient(low = "blue", high = "red") +  # couleur = moyenne d'expression
  labs(
    x = "Genes",
    y = "Patterns",
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

# Ajouter un contour blanc sans changer la taille des ronds
p$layers[[1]]$aes_params$stroke <- 0.5  # Contour blanc autour des ronds
p$layers[[1]]$aes_params$color <- "white"  # Contour blanc
p$layers[[1]]$aes_params$fill <- p$layers[[1]]$aes_params$fill  # Laisser la couleur de remplissage inchangée

# Afficher
print(p)


#####Meme plot mais avec Vim
# Gènes à visualiser
features <- c("SV40", "TdTomato", "Krt14", "Vim")

# Générer le DotPlot standard
p <- DotPlot(seurat_obj, features = features) + RotatedAxis()

# Modifier uniquement le thème et les contours sans changer les échelles
p <- p +
  scale_color_gradient(low = "blue", high = "red") +  # couleur = moyenne d'expression
  labs(
    x = "Genes",
    y = "Patterns",
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

# Ajouter un contour blanc sans changer la taille des ronds
p$layers[[1]]$aes_params$stroke <- 0.5  # Contour blanc autour des ronds
p$layers[[1]]$aes_params$color <- "white"  # Contour blanc
p$layers[[1]]$aes_params$fill <- p$layers[[1]]$aes_params$fill  # Laisser la couleur de remplissage inchangée

# Afficher
print(p)
