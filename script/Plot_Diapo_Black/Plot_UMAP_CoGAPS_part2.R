# --- Librairies nécessaires
library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)

# --- Chargement des données
seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")
DefaultAssay(seurat_obj) <- "RNA" 
seurat_obj$pattern_cogaps <- gsub("Pattern_", "", seurat_obj$pattern_cogaps)
seurat_obj$pattern_cogaps <- factor(seurat_obj$pattern_cogaps)  # garder ordre original

# --- Charger la signature ImSig
signature_imsig <- read.table("C:/Users/louis/Desktop/Stage/scrnaseq/signatures/ImSig2_mouse.txt", header = FALSE, stringsAsFactors = FALSE)$V1

# --- Calcul du score de module ImSig
seurat_obj <- AddModuleScore(seurat_obj, features = list(signature_imsig), name = "ImSig_Score")

# --- Définir palette de couleurs UMAP pour clusters (couleurs inversées uniquement)
cluster_levels <- levels(seurat_obj$pattern_cogaps)
umap_colors <- scales::hue_pal()(length(cluster_levels))
umap_colors <- rev(umap_colors)  # inverser l'ordre des couleurs seulement
names(umap_colors) <- cluster_levels

# --- DotPlot : TdTomato, Krt14, Vim
genes_of_interest <- c("TdTomato", "Krt14", "Vim")

pC <- DotPlot(seurat_obj, features = genes_of_interest, group.by = "pattern_cogaps") +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_family = "sans") +
  labs(
    x = "Genes",
    y = "CoGAPS patterns"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "white", size = 14, face = "bold"),
    axis.text.y = element_text(
      size = 16,
      face = "bold",
      color = umap_colors[cluster_levels]  # couleurs inversées mais labels dans l'ordre
    ),
    axis.title.x = element_text(color = "white", size = 16, margin = margin(t = 5)),
    axis.title.y = element_text(color = "white", size = 16),
    axis.line = element_line(color = "white"),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    plot.title = element_blank(),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

# --- Fonction pour annoter le nombre de cellules
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
    size = 3,
    fontface = "bold",
    inherit.aes = FALSE
  )
}

# --- Violin Plot du score ImSig
pD <- VlnPlot(
  object = seurat_obj,
  features = "ImSig_Score1",
  group.by = "pattern_cogaps",
  pt.size = 0
) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
  labs(y = "ImSig Score", x = "CoGAPS Patterns") +
  theme_minimal(base_family = "sans") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, color = "white"),
    axis.text.y = element_text(size = 16, color = "white"),
    axis.title.x = element_text(color = "white", size = 16, margin = margin(t = 5)),
    axis.title.y = element_text(color = "white", size = 16),
    axis.line = element_line(color = "white"),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    plot.title = element_blank()
  ) +
  scale_fill_manual(values = umap_colors)  # <-- appliquer la palette inversée ici

# Ajustement de l'échelle Y
max_y1 <- max(seurat_obj@meta.data$ImSig_Score1[!is.na(seurat_obj@meta.data$pattern_cogaps)], na.rm = TRUE)
pD <- pD + scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(NA, max_y1 * 1.15))

# Ajout des annotations de cellules
pD <- add_cell_count(pD, seurat_obj@meta.data, "pattern_cogaps")

# --- UMAP avec labels blancs et plus gros
umap <- DimPlot(seurat_obj, reduction = "umap", group.by = "pattern_cogaps", label = TRUE, pt.size = 0.5, label.size = 8, label.color = "white") +
  scale_color_manual(values = umap_colors, name = "Clusters") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(title = NULL) +
  theme(
    plot.title = element_blank(),
    legend.text = element_text(size = 12, color = "white"),
    axis.text = element_text(size = 13, color = "white"),
    axis.title = element_text(size = 15, color = "white"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white")
  )

# --- Assemblage des deux plots (DotPlot + Violin) avec tags C et D
panel_plot <- pC + pD + 
  plot_layout(ncol = 2, widths = c(1, 1)) + 
  theme(
    plot.background = element_rect(fill = "black", color = NA)
  )

# --- Sauvegarde du panel DotPlot + Violin
ggsave(
  "C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Diapo_Black/Panel_ImSig_DotPlot_CoGAPS.png",
  plot = panel_plot, width = 13, height = 6, dpi = 300, bg = "black"
)

