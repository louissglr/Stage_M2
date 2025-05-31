library(Seurat)
library(dplyr)
library(ggplot2)
library(ggtext)  # Pour element_markdown()

seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")
seurat_obj <- SetIdent(seurat_obj, value = seurat_obj$pattern_cogaps)
DefaultAssay(seurat_obj) <- "RNA"

pattern_order <- c(
  "Pattern_1", "Pattern_4", "Pattern_5",        # Immune
  "Pattern_6",                                  # Pulmonary
  "Pattern_2", "Pattern_3", "Pattern_7", "Pattern_8", "Pattern_9"  # Cancerous
)
seurat_obj@active.ident <- factor(seurat_obj@active.ident, levels = pattern_order)

markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

non_coding_pattern <- "^Gm|Rik|Mir|^mt-|^Rp[sl]|^Snr|^Malat1|^Xist"
markers_filtered <- markers %>%
  filter(!grepl(non_coding_pattern, gene, ignore.case = TRUE))

top3 <- markers_filtered %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC)

genes_top <- unique(top3$gene)
genes_to_plot <- c(setdiff(genes_top, c("TdTomato", "Krt14")), "Krt14", "TdTomato")

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

# Créer les labels Y colorés en HTML
colored_labels <- sapply(levels(seurat_obj@active.ident), function(x) {
  color <- pattern_colors[x]
  paste0("<span style='color:", color, "'>", x, "</span>")
})
names(colored_labels) <- levels(seurat_obj@active.ident)

# Créer les labels X colorés en HTML (Krt14 et TdTomato en rouge)
labels_x <- sapply(genes_to_plot, function(gene) {
  if (gene %in% c("Krt14", "TdTomato")) {
    paste0("<span style='color:red; font-weight:bold'>", gene, "</span>")
  } else {
    gene
  }
})

# DotPlot
p <- DotPlot(seurat_obj, features = genes_to_plot) +
  RotatedAxis() +
  scale_x_discrete(labels = labels_x) +  # appliquer labels x colorés
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1, size = 11, face = "bold"),  # texte x en markdown pour interpréter HTML
    axis.text.y = element_markdown(size = 11, face = "bold"),  # texte y en markdown
    plot.margin = margin(10, 10, 10, 50)
  ) +
  scale_y_discrete(labels = colored_labels)

# Ajouter la bande des catégories à gauche (optionnel, à garder si tu veux)
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
            size = ifelse(annot_data$category == "Pulmonary", 3, 4.5),
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

ggsave(filename = "C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Rapport/dotplot_cogaps.png", 
       plot = final_plot, width = 10, height = 8, units = "in", dpi = 300, bg = "white")
