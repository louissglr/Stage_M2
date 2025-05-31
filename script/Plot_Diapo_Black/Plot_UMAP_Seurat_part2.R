library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)

# --- Chargement des données
seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds")
signature_imsig <- read.table("C:/Users/louis/Desktop/Stage/scrnaseq/signatures/ImSig2_mouse.txt", header = FALSE, stringsAsFactors = FALSE)$V1

# --- Calcul du score de module ImSig
seurat_obj <- AddModuleScore(seurat_obj, features = list(signature_imsig), name = "ImSig_Score")

# --- Définir palette de couleurs UMAP pour clusters
cluster_levels <- levels(seurat_obj$seurat_clusters)
umap_colors <- scales::hue_pal()(length(cluster_levels))
names(umap_colors) <- cluster_levels

# --- A. DotPlot : TdTomato, Krt14, Vim
genes_of_interest <- c("TdTomato", "Krt14", "Vim")

pA <- DotPlot(seurat_obj, features = genes_of_interest, group.by = "seurat_clusters") +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 16) +
  labs(
    x = "Genes",
    y = "Seurat clusters"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "white", size = 14),
    axis.text.y = element_text(
      size = 16,
      face = "bold",
      color = umap_colors[cluster_levels]
    ),
    axis.title.x = element_text(color = "white", size = 16, margin = margin(t = 5)),
    axis.title.y = element_text(color = "white", size = 16),
    axis.line = element_line(color = "white"),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    
    legend.background = element_rect(fill = NA, color = NA),
    legend.key = element_rect(fill = NA, color = NA),
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
    size = 4,
    fontface = "bold",
    inherit.aes = FALSE
  )
}

# --- B. Violin Plot du score ImSig
pB <- VlnPlot(
  object = seurat_obj,
  features = "ImSig_Score1",
  group.by = "seurat_clusters",
  pt.size = 0
) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
  labs(y = "ImSig Score", x = "Seurat clusters") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, color = "white"),
    axis.text.y = element_text(size = 16, color = "white"),
    axis.title.x = element_text(color = "white", size = 16, margin = margin(t = 5)),
    axis.title.y = element_text(color = "white", size = 16),
    axis.line = element_line(color = "white"),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA)
  )

# Ajustement de l'échelle Y
max_y1 <- max(seurat_obj@meta.data$ImSig_Score1[!is.na(seurat_obj@meta.data$seurat_clusters)], na.rm = TRUE)
pB <- pB + scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(NA, max_y1 * 1.15))

# Ajout des annotations de cellules
pB <- add_cell_count(pB, seurat_obj@meta.data, "seurat_clusters")

# --- Assemblage des deux plots (panel A + B) sans tags
panel_plot <- pA + pB + 
  plot_layout(ncol = 2, widths = c(1, 1)) + 
  theme(
    plot.background = element_rect(fill = "black", color = NA)
  )

# --- Sauvegarde du panel
ggsave("C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Rapport/Panel_ImSig_DotPlot_noTags.png",
       plot = panel_plot, width = 12, height = 6, dpi = 300, bg = "black")


# --- Sauvegarde du panel
ggsave("C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Rapport/Panel_ImSig_DotPlot_black.png",
       plot = panel_plot, width = 12, height = 6, dpi = 300, bg = "black")


# --- Sauvegarde du panel
ggsave("C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Diapo_Black/Panel_ImSig_DotPlot_Seurat.png",
       plot = panel_plot, width = 12, height = 6, dpi = 300, bg = "black")
