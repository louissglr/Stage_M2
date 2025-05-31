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

# --- C. DotPlot : TdTomato, Krt14, Vim
genes_of_interest <- c("TdTomato", "Krt14")

pC <- DotPlot(seurat_obj, features = genes_of_interest, group.by = "pattern_cogaps") +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    x = "Genes",
    y = "CoGAPS patterns"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
    axis.text.y = element_text(
      size = 16,
      face = "bold",
      color = umap_colors[cluster_levels]  # couleurs inversées mais labels dans l'ordre
    ),
    axis.title.x = element_text(color = "black", size = 16, margin = margin(t = 5)),
    axis.title.y = element_text(color = "black", size = 16),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_blank()
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
    color = "black",
    size = 2.7,
    fontface = "bold",
    inherit.aes = FALSE
  )
}

# --- D. Violin Plot du score ImSig
pD <- VlnPlot(
  object = seurat_obj,
  features = "ImSig_Score1",
  group.by = "pattern_cogaps",
  pt.size = 0
) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
  labs(y = "ImSig Score", x = "CoGAPS Patterns") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(color = "black", size = 16, margin = margin(t = 5)),
    axis.title.y = element_text(color = "black", size = 16),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_blank()
  ) +
  scale_fill_manual(values = umap_colors)  # <-- appliquer la palette inversée ici

# Ajustement de l'échelle Y
max_y1 <- max(seurat_obj@meta.data$ImSig_Score1[!is.na(seurat_obj@meta.data$pattern_cogaps)], na.rm = TRUE)
pD <- pD + scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(NA, max_y1 * 1.15))

# Ajout des annotations de cellules
pD <- add_cell_count(pD, seurat_obj@meta.data, "pattern_cogaps")

# --- Assemblage des deux plots (panel C + D) avec tags C et D
panel_plot <- pC + pD + 
  plot_layout(ncol = 2, widths = c(1, 1)) + 
  plot_annotation(tag_levels = list(c("C", "D"))) & 
  theme(plot.tag = element_text(face = "bold", size = 20, hjust = 0, vjust = 1))

# --- Sauvegarde du panel
ggsave("C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Rapport/Panel_ImSig_DotPlot_CoGAPS.png",
       plot = panel_plot, width = 13, height = 6, dpi = 300)
