# Charger les bibliothèques nécessaires
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)  # Pour assembler les plots

# Fonctions utiles ---------------------------------------------------------

# Préparation des données avec possibilité de forcer min/max pour certaines variables
prepare_data <- function(repertoire, var_order, ident_order, global_minmax = NULL) {
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
  
  # df_blank pour percent.ribo et percent.mt forcés entre 0 et 100
  df_blank <- df_long %>%
    filter(variable %in% c("percent.ribo", "percent.mt")) %>%
    group_by(orig.ident, variable) %>%
    summarise(min_val = 0, max_val = 100, .groups = "drop")
  
  # Ajouter les min/max globaux pour nCount_RNA et nFeature_RNA si fournis
  if (!is.null(global_minmax)) {
    df_blank <- bind_rows(df_blank, global_minmax)
  }
  
  return(list(df_long = df_long, seuils = seuils, df_blank = df_blank))
}

# Fonction pour tracer le graphique QC clair (sans thème noir)
plot_QC <- function(df_long, seuils, df_blank) {
  ggplot(df_long, aes(x = orig.ident, y = value, fill = orig.ident)) +
    geom_boxplot(outlier.shape = NA, coef = 10, color = "black") +
    geom_blank(data = df_blank, aes(y = min_val)) +   # Forcer minimum
    geom_blank(data = df_blank, aes(y = max_val)) +   # Forcer maximum
    theme_minimal(base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      strip.text = element_text(hjust = 0.5, color = "black", size = 12, face = "bold"),
      axis.title.x = element_text(color = "black", hjust = 0.5, size = 17),
      axis.title.y = element_text(color = "black", hjust = 0.5, size = 17),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(hjust = 0.5, color = "black", size = 16, face = "bold")
    ) +
    labs(x = NULL, y = NULL) +
    scale_fill_brewer(name = "Samples", palette = "Set3") +
    facet_wrap(vars(variable), scales = "free_y") +
    geom_hline(data = seuils, aes(yintercept = seuil, color = color), linetype = "dashed", size = 0.5) +
    scale_color_manual(name = "Threshold QC", values = c("red" = "red"))
}

# Calcul min/max global ----------------------------------------------------

var_order <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "doublets.score", "log10GenesPerUMI")
ident_order <- c("N1-0h", "N1-72h", "N1-1week", "N1-2week", "N1-4week", "N2-0h", "N2-72h", "N2-1week")

get_minmax <- function(repertoire, var) {
  fichiers_rds <- list.files(repertoire, pattern = "\\.rds$", full.names = TRUE)
  objets_seurat <- lapply(fichiers_rds, readRDS)
  seurat_object <- merge(objets_seurat[[1]], y = objets_seurat[-1], merge.data = TRUE)
  seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
  meta <- seurat_object@meta.data
  res <- meta %>% summarise(min_val = min(!!sym(var)), max_val = max(!!sym(var)))
  return(res)
}

# Pour nCount_RNA
minmax_nCount_final <- get_minmax("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/02_NormalizedData", "nCount_RNA")
minmax_nCount_full <- get_minmax("C:/Users/louis/Desktop/Stage/scrnaseq/output/02_NormalizedData", "nCount_RNA")

min_nCount <- min(minmax_nCount_final$min_val, minmax_nCount_full$min_val)
max_nCount <-  max(minmax_nCount_final$max_val, minmax_nCount_full$max_val)

# Pour nFeature_RNA
minmax_nFeature_final <- get_minmax("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/02_NormalizedData", "nFeature_RNA")
minmax_nFeature_full <- get_minmax("C:/Users/louis/Desktop/Stage/scrnaseq/output/02_NormalizedData", "nFeature_RNA")

min_nFeature <- min(minmax_nFeature_final$min_val, minmax_nFeature_full$min_val)
max_nFeature <- max(minmax_nFeature_final$max_val, minmax_nFeature_full$max_val)

# Dataframe pour forcer min/max dans geom_blank
global_minmax <- data.frame(
  orig.ident = factor(rep(ident_order, 2), levels = ident_order),
  variable = factor(rep(c("nCount_RNA", "nFeature_RNA"), each = length(ident_order)), levels = var_order),
  min_val = c(rep(min_nCount, length(ident_order)), rep(min_nFeature, length(ident_order))),
  max_val = c(rep(max_nCount, length(ident_order)), rep(max_nFeature, length(ident_order)))
)

# Préparer les données -----------------------------------------------------

data_gene <- prepare_data("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/02_NormalizedData", var_order, ident_order, global_minmax)
data_fullgene <- prepare_data("C:/Users/louis/Desktop/Stage/scrnaseq/output/02_NormalizedData", var_order, ident_order, global_minmax)

# Créer les graphiques ------------------------------------------------------

p_gene <- plot_QC(data_gene$df_long, data_gene$seuils, data_gene$df_blank)
p_fullgene <- plot_QC(data_fullgene$df_long, data_fullgene$seuils, data_fullgene$df_blank)

# Assembler en figure finale ------------------------------------------------
figure_finale <- p_gene / p_fullgene + 
  plot_layout(heights = c(1,1), guides = "collect") +  # <- fusion des légendes
  plot_annotation(
    tag_levels = list(c("A", "B")),
    tag_prefix = "",
    tag_suffix = ""
  ) & 
  theme(
  plot.tag = element_text(size = 24, face = "bold", family = "Arial"),
  plot.tag.position = c(0.02, 0.98),
  legend.position = "bottom"
)



# Sauvegarder la figure ----------------------------------------------------
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/03_FilteredData/quality_controls_plot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(file.path(output_dir, "figure_QC_Gene_FullGene_same_scale.png"), plot = figure_finale, width = 12, height = 14, dpi = 300)

# Afficher la figure --------------------------------------------------------
print(figure_finale)

