# Charger l'objet intégré
int <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/15_IntegrationWithPerou/integration.rds")

# Créer le répertoire de sortie s’il n’existe pas
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/16_Analyze_Integration"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

library(Seurat)
library(ggplot2)

# Liste des réductions UMAP à afficher
reductions <- c("umap.cca", "umap.rpca", "umap.harmony")

# Boucle avec thème fond noir
for (reduc in reductions) {
  p <- DimPlot(int, reduction = reduc, group.by = "seurat_clusters", pt.size = 0.3) +
    ggtitle(paste("UMAP -", reduc, "- orig.ident")) +
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
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      plot.title = element_text(color = "white", hjust = 0.5)
    ) +
    coord_fixed()
  
  ggsave(
    filename = file.path(output_dir, paste0("umap_", reduc, "_orig_ident_black.png")),
    plot = p, bg = "black", width = 7, height = 6, dpi = 300
  )
  print(p)
}


# Annoter les cellules
int$cell_type_simple <- ifelse(
  grepl("^N1|^N2", Cells(int)),
  "exp",
  "ref"
)
int$cell_type_simple <- factor(int$cell_type_simple, levels = c("ref", "exp"))

# Fonction de plot avec thème noir
plot_umap_black <- function(obj, reduction_name, output_name) {
  p <- DimPlot(
    object = obj,
    reduction = reduction_name,
    group.by = "cell_type_simple",
    cols = c("ref" = "white", "exp" = "red"),
    pt.size = 0.2
  ) +
    ggtitle(paste("UMAP -", reduction_name)) +
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
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      plot.title = element_text(color = "white", hjust = 0.5)
    ) +
    coord_fixed()
  
  ggsave(file.path(output_dir, output_name), plot = p, bg = "black", width = 7, height = 6, dpi = 300)
  print(p)
}

# Générer les UMAPs noirs
plot_umap_black(int, "umap.cca", "umap_cca_ref_vs_exp_black.png")
plot_umap_black(int, "umap.rpca", "umap_rpca_ref_vs_exp_black.png")
plot_umap_black(int, "umap.harmony", "umap_harmony_ref_vs_exp_black.png")

library(scales)

# Annoter les groupes personnalisés
int$cell_group_custom <- "other"
int$cell_group_custom[int$seurat_clusters %in% c(2, 3, 9)] <- "cluster_group"
int$cell_group_custom[int$md_cluster == "immune"] <- "immune_group"

# Palette personnalisée
color_map_custom <- c(
  "cluster_group" = alpha("blue", 0.8),
  "immune_group" = alpha("green", 0.8),
  "other" = alpha("white", 0.2)
)

# Fonction pour UMAP custom
plot_umap_custom <- function(obj, reduction_name, output_name) {
  p <- DimPlot(
    object = obj,
    reduction = reduction_name,
    group.by = "cell_group_custom",
    cols = color_map_custom,
    pt.size = 0.3
  ) +
    ggtitle(paste("UMAP -", reduction_name)) +
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
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      plot.title = element_text(color = "white", hjust = 0.5)
    ) +
    coord_fixed()
  
  ggsave(file.path(output_dir, output_name), plot = p, bg = "black", width = 7, height = 6, dpi = 300)
  print(p)
}

# UMAPs personnalisés
plot_umap_custom(int, "umap.cca", "custom_umap_cca.png")
plot_umap_custom(int, "umap.rpca", "custom_umap_rpca.png")
plot_umap_custom(int, "umap.harmony", "custom_umap_harmony.png")

# Ajouter les scores ImSig
signature_imsig <- read.table("C:/Users/louis/Desktop/Stage/scrnaseq/signatures/ImSig2_mouse.txt", header = FALSE, stringsAsFactors = FALSE)$V1
int <- JoinLayers(int)
int <- AddModuleScore(int, features = list(signature_imsig), name = "ImSig_Score")

# Gradient rouge
gradient_colors <- c("white", "red")

# Fonction pour FeaturePlot stylé
plot_feature_custom <- function(obj, feature_name, reduction_name, output_name) {
  p <- FeaturePlot(
    object = obj,
    features = feature_name,
    reduction = reduction_name,
    pt.size = 0.3,
    cols = gradient_colors
  ) +
    ggtitle(paste("UMAP -", reduction_name, "-", feature_name)) +
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
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      plot.title = element_text(color = "white", hjust = 0.5)
    ) +
    coord_fixed()
  
  ggsave(file.path(output_dir, output_name), plot = p, bg = "black", width = 7, height = 6, dpi = 300)
  print(p)
}

# Générer les plots ImSig
plot_feature_custom(int, "ImSig_Score1", "umap.cca", "umap_imsig_cca.png")
plot_feature_custom(int, "ImSig_Score1", "umap.rpca", "umap_imsig_rpca.png")
plot_feature_custom(int, "ImSig_Score1", "umap.harmony", "umap_imsig_harmony.png")
