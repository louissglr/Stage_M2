# Charger les bibliothèques nécessaires
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggh4x)

# Définir l'ordre des échantillons et la palette commune ------------------

ordre_legende <- c("N1-0h", "N1-72h", "N1-1week", "N1-2week", "N1-4week", 
                   "N2-0h", "N2-72h", "N2-1week")

palette_orig <- setNames(brewer.pal(n = length(ordre_legende), "Set2"), ordre_legende)

# Fonctions utiles ---------------------------------------------------------

prepare_data <- function(repertoire, var_order, ident_order) {
  fichiers_rds <- list.files(repertoire, pattern = "\\.rds$", full.names = TRUE)
  objets_seurat <- lapply(fichiers_rds, readRDS)
  seurat_object <- merge(objets_seurat[[1]], y = objets_seurat[-1], merge.data = TRUE)
  seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
  
  seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)
  
  df_long <- seurat_object@meta.data %>%
    select(orig.ident, all_of(var_order)) %>%
    pivot_longer(cols = all_of(var_order), names_to = "variable", values_to = "value") %>%
    mutate(
      variable = factor(variable, levels = var_order),
      orig.ident = factor(orig.ident, levels = ident_order)
    ) %>%
    arrange(variable)
  
  seuils <- data.frame(
    variable = factor(
      c("nCount_RNA", "nCount_RNA", "nFeature_RNA", "nFeature_RNA", 
        "percent.mt", "doublets.score", "log10GenesPerUMI"),
      levels = var_order
    ),
    seuil = c(50, 7000, 200, 5000, 15, 0.6, 0.7),
    color = "red"
  )
  
  # df_blank avec min/max manuels
  df_blank <- df_long %>%
    filter(variable %in% c("percent.ribo", "percent.mt")) %>%
    group_by(orig.ident, variable) %>%
    summarise(min_val = 0, max_val = 100, .groups = "drop")
  
  # Ajouter min/max manuels pour nCount_RNA et nFeature_RNA
  df_blank <- bind_rows(
    df_blank,
    expand.grid(orig.ident = ident_order, variable = "nCount_RNA") %>%
      mutate(min_val = 0, max_val = 10000),
    expand.grid(orig.ident = ident_order, variable = "nFeature_RNA") %>%
      mutate(min_val = 0, max_val = NA)  # facultatif ici
  )
  
  return(list(df_long = df_long, seuils = seuils, df_blank = df_blank))
}

# Fonction de tracé avec ggh4x pour fixer les échelles -----------------------

plot_QC <- function(df_long, seuils, df_blank, palette_orig) {
  ggplot(df_long, aes(x = orig.ident, y = value, fill = orig.ident)) +
    geom_boxplot(outlier.shape = NA, coef = 10, color = "black") +
    geom_blank(data = df_blank, aes(y = min_val)) +
    geom_blank(data = df_blank, aes(y = max_val)) +
    theme_minimal(base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.ticks = element_line(color = "black"),
      strip.text = element_text(hjust = 0.5, color = "black", size = 13, face = "bold"),
      axis.title.x = element_text(color = "black", hjust = 0.5, size = 18),
      axis.title.y = element_text(color = "black", hjust = 0.5, size = 18),
      legend.text = element_text(color = "black", size = 12),
      legend.title = element_text(color = "black", size = 13),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(hjust = 0.5, color = "black", size = 17, face = "bold")
    ) +
    labs(x = NULL, y = NULL) +
    scale_fill_manual(name = "Samples", values = palette_orig) +
    facet_wrap(vars(variable), scales = "free_y") +
    geom_hline(data = seuils, aes(yintercept = seuil, color = color), linetype = "dashed", size = 0.5) +
    scale_color_manual(name = "Threshold QC", values = c("red" = "red")) +
    ggh4x::facetted_pos_scales(
      y = list(
        `nCount_RNA` = scale_y_continuous(limits = c(0, 10000))
      )
    )
}

# Définition des variables et ordre des échantillons -----------------------

var_order <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "doublets.score", "log10GenesPerUMI")
ident_order <- ordre_legende

# Préparer les données -----------------------------------------------------

data_gene <- prepare_data("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/02_NormalizedData", var_order, ident_order)

# Créer le graphique -------------------------------------------------------

p_gene <- plot_QC(data_gene$df_long, data_gene$seuils, data_gene$df_blank, palette_orig)

# Sauvegarder le graphique -------------------------------------------------

output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/03_FilteredData/quality_controls_plot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(
  file.path(output_dir, "figure_QC_Gene_only.png"),
  plot = p_gene,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Afficher le graphique ----------------------------------------------------
print(p_gene)

# Calcul des médianes ------------------------------------------------------

get_medians <- function(df_long) {
  df_long %>%
    group_by(orig.ident, variable) %>%
    summarise(mediane = median(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(variable, orig.ident)
}

medians_gene <- get_medians(data_gene$df_long)

print("Médianes pour les données GENE :")
print(medians_gene)

write.csv(medians_gene, file = file.path(output_dir, "medians_gene.csv"), row.names = FALSE)
