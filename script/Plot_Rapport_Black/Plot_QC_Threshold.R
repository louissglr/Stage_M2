# Charger les bibliothèques
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)

# Charger les fichiers Seurat
repertoire <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/02_NormalizedData"
fichiers_rds <- list.files(repertoire, pattern = "\\.rds$", full.names = TRUE)
objets_seurat <- lapply(fichiers_rds, readRDS)

# Fusionner les objets Seurat
seurat_object <- merge(objets_seurat[[1]], y = objets_seurat[-1], merge.data = TRUE)
seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])

# Calculer log10GenesPerUMI
seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)

# Définir les ordres
var_order <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "doublets.score", "log10GenesPerUMI")
ident_order <- c("N1-0h", "N1-72h", "N1-1week", "N1-2week", "N1-4week", "N2-0h", "N2-72h", "N2-1week")

# Mettre les métadonnées au format long
df_long <- seurat_object@meta.data %>%
  select(orig.ident, all_of(var_order)) %>%
  pivot_longer(cols = all_of(var_order), names_to = "variable", values_to = "value") %>%
  mutate(
    variable = factor(variable, levels = var_order),
    orig.ident = factor(orig.ident, levels = ident_order)
  ) %>%
  arrange(variable)

# Définir les seuils QC à afficher
seuils <- data.frame(
  variable = factor(
    c("nCount_RNA", "nCount_RNA", "nFeature_RNA", "nFeature_RNA", 
      "percent.mt", "doublets.score", "log10GenesPerUMI"),
    levels = var_order
  ),
  seuil = c(50, 7000, 200, 5000, 15, 0.6, 0.7),
  color = "red"
)

# Créer le dataframe pour forcer l’échelle 0-100 pour percent.ribo et percent.mt
df_blank <- df_long %>%
  filter(variable %in% c("percent.ribo", "percent.mt")) %>%
  group_by(orig.ident, variable) %>%
  summarise(min_val = 0, max_val = 100, .groups = "drop")

# Créer le graphique QC
p <- ggplot(df_long, aes(x = orig.ident, y = value, fill = orig.ident)) +
  geom_boxplot(outlier.shape = NA, coef = 10, color = "white") +
  geom_blank(data = df_blank, aes(y = min_val)) +   # Forcer minimum 0
  geom_blank(data = df_blank, aes(y = max_val)) +   # Forcer maximum 100
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "white"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "white"),
    strip.text = element_text(hjust = 0.5, color = "white"),
    axis.title.x = element_text(color = "white", hjust = 0.5, size = 16),
    axis.title.y = element_text(color = "white", hjust = 0.5, size = 16),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "transparent"),
    panel.grid.minor = element_line(color = "transparent"),
    panel.border = element_rect(color = "white", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, color = "white", size = 18, face = "bold")
  ) +
  labs(x = NULL, y = NULL) +
  scale_fill_brewer(name = "Samples", palette = "Set3") +
  facet_wrap(vars(variable), scales = "free_y") +
  geom_hline(data = seuils, aes(yintercept = seuil, color = color), linetype = "dashed", size = 0.5) +
  scale_color_manual(name = "Threshold QC", values = c("red" = "red"))

# Créer le répertoire de sortie s’il n’existe pas
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/03_FilteredData/quality_controls_plot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Sauvegarder le graphique
ggsave(filename = file.path(output_dir, "graphique_QC_complet.png"), plot = p, width = 10, height = 8, dpi = 300)

# Afficher le graphique
print(p)


# Appliquer tous les filtres QC combinés
seurat_qc_filtree <- subset(
  seurat_object,
  subset = 
    nCount_RNA >= 50 &
    nCount_RNA <= 7000 &
    nFeature_RNA >= 200 &
    nFeature_RNA <= 5000 &
    percent.mt <= 15 &
    doublets.score <= 0.6 &
    log10GenesPerUMI >= 0.7
)

# Compter les cellules
total_cellules <- ncol(seurat_object)
cellules_restantes <- ncol(seurat_qc_filtree)
cellules_filtrees <- total_cellules - cellules_restantes

# Résumé dans un dataframe
df_resume_global <- data.frame(
  total_cellules = total_cellules,
  cellules_restantes = cellules_restantes,
  cellules_filtrees = cellules_filtrees
)

# Afficher
print(df_resume_global)


