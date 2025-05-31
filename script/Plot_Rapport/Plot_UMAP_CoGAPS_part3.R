library(Seurat)
library(dplyr)
library(ggplot2)
library(ggtext)      # pour element_markdown()
library(patchwork)

# --- Fonction pour créer le DotPlot amélioré (final_plot)
create_final_plot <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- SetIdent(seurat_obj, value = seurat_obj$pattern_cogaps)
  
  # Ordre des patterns
  pattern_order <- c(
    "Pattern_1", "Pattern_4", "Pattern_5",        # Immune
    "Pattern_6",                                  # Pulmonary
    "Pattern_2", "Pattern_3", "Pattern_7", "Pattern_8", "Pattern_9"  # Cancerous
  )
  seurat_obj@active.ident <- factor(seurat_obj@active.ident, levels = pattern_order)
  
  # Trouver marqueurs (top 3 par cluster) — adapte si tu veux refaire cette étape
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  non_coding_pattern <- "^Gm|Rik|Mir|^mt-|^Rp[sl]|^Snr|^Malat1|^Xist"
  markers_filtered <- markers %>%
    filter(!grepl(non_coding_pattern, gene, ignore.case = TRUE))
  top3 <- markers_filtered %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC)
  
  genes_top <- unique(top3$gene)
  genes_to_plot <- c(setdiff(genes_top, c("TdTomato", "Krt14", "Fos")), "Krt14", "TdTomato")
  
  
  # Couleurs des patterns
  pattern_colors <- c(
    "Pattern_9" = "#f8776c",
    "Pattern_8" = "#d39301",
    "Pattern_7" = "#92ab00",
    "Pattern_6" = "#00ba38",
    "Pattern_5" = "#01c09f",
    "Pattern_4" = "#00b8e2",
    "Pattern_3" = "#72a6ff",
    "Pattern_2" = "#db73fa",
    "Pattern_1" = "#ff60c2"
  )
  
  # Labels Y colorés en HTML
  colored_labels <- sapply(levels(seurat_obj@active.ident), function(x) {
    color <- pattern_colors[x]
    paste0("<span style='color:", color, "'>", x, "</span>")
  })
  names(colored_labels) <- levels(seurat_obj@active.ident)
  
  # Labels X colorés en HTML (Krt14 et TdTomato en rouge)
  labels_x <- sapply(genes_to_plot, function(gene) {
    if (gene %in% c("Krt14", "TdTomato")) {
      paste0("<span style='color:red; font-weight:bold'>", gene, "</span>")
    } else {
      gene
    }
  })
  
  # DotPlot de base avec échelle de couleur continue bleu->rouge
  p <- DotPlot(seurat_obj, features = genes_to_plot) +
    RotatedAxis() +
    scale_x_discrete(labels = labels_x) +
    scale_y_discrete(labels = colored_labels) +
    scale_color_gradient(low = "blue", high = "red") +   # <- Palette couleur continue
    theme(
      axis.text.x = element_markdown(angle = 45, hjust = 1, size = 7, face = "bold"),
      axis.text.y = element_markdown(size = 14, face = "bold"),
      plot.margin = margin(10, 10, 10, 50)
    )
  
  # Bande des catégories à gauche
  annot_data <- data.frame(
    category = c("Immune", "Pulmonary", "Cancerous"),
    ymin = c(0.5, 3.5, 4.5),
    ymax = c(3.5, 4.5, 9.5),
    color = c("lightblue", "lightgreen", "lightpink")
  )
  
  final_plot <- p +
    geom_rect(data = annot_data, aes(xmin = -0.5, xmax = 0, ymin = ymin, ymax = ymax, fill = category),
              inherit.aes = FALSE, alpha = 0.3) +
    geom_text(data = annot_data, aes(x = -0.25, y = (ymin + ymax)/2, label = category),
              inherit.aes = FALSE,
              size = ifelse(annot_data$category == "Pulmonary", 2.5, 3),
              fontface = "bold", angle = 90) +
    scale_fill_manual(
      values = c(
        "Immune" = "lightblue",
        "Pulmonary" = "lightgreen",
        "Cancerous" = "lightpink"
      ),
      guide = "none"
    ) +
    coord_cartesian(xlim = c(-1, length(unique(p$data$features.plot))))
  
  return(final_plot)
}

# --- Chargement des données
seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")
DefaultAssay(seurat_obj) <- "RNA" 
seurat_obj$pattern_cogaps <- gsub("Pattern_", "Pattern_", seurat_obj$pattern_cogaps)  # Remettre le prefix si besoin
seurat_obj$pattern_cogaps <- factor(seurat_obj$pattern_cogaps)  # garder ordre original

# --- Créer final_plot (DotPlot coloré + annotations)
final_plot <- create_final_plot(seurat_obj)

# --- Charger signature ImSig et calcul du score
signature_imsig <- read.table("C:/Users/louis/Desktop/Stage/scrnaseq/signatures/ImSig2_mouse.txt", header = FALSE, stringsAsFactors = FALSE)$V1
seurat_obj <- AddModuleScore(seurat_obj, features = list(signature_imsig), name = "ImSig_Score")

# --- Palette inversée pour UMAP / patterns
cluster_levels <- levels(seurat_obj$pattern_cogaps)
umap_colors <- scales::hue_pal()(length(cluster_levels))
umap_colors <- rev(umap_colors)
names(umap_colors) <- cluster_levels

# --- Violin Plot du score ImSig
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
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(color = "black", size = 16, margin = margin(t = 5)),
    axis.title.y = element_text(color = "black", size = 16),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_blank()
  ) +
  scale_fill_manual(values = umap_colors)

# Ajuster échelle Y
max_y1 <- max(seurat_obj@meta.data$ImSig_Score1[!is.na(seurat_obj@meta.data$pattern_cogaps)], na.rm = TRUE)
pD <- pD + scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(NA, max_y1 * 1.15))

# Fonction pour ajouter le nombre de cellules
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

pD <- add_cell_count(pD, seurat_obj@meta.data, "pattern_cogaps")

# --- Assemblage des plots : final_plot (DotPlot amélioré) + pD (Violin)
panel_plot <- final_plot + pD + 
  plot_layout(ncol = 2, widths = c(1, 1)) + 
  plot_annotation(tag_levels = list(c("C", "D"))) & 
  theme(plot.tag = element_text(face = "bold", size = 20, hjust = 0, vjust = 1))

# --- Sauvegarde
ggsave("C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Rapport/Panel_ImSig_DotPlot_CoGAPS2.png",
       plot = panel_plot, width = 13, height = 6, dpi = 300)
