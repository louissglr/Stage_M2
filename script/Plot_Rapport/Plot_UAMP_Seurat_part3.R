# Chargement des packages
library(dplyr)
library(ggplot2)
library(ggtext)
library(patchwork)
library(Seurat)

# Chargement de l'objet Seurat
seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")
seurat_obj <- SetIdent(seurat_obj, value = seurat_obj$seurat_clusters)
DefaultAssay(seurat_obj) <- "RNA"

# Ordre des clusters avec 11 dans Unknown, 12 dans Immune
pattern_order <- c(
  "2", "3", "4", "7", "9", "10", "12",    # Immune
  "6", "8",                               # Pulmonary
  "0", "1", "5",                          # Cancerous
  "11", "13"                              # Unknown
)
seurat_obj@active.ident <- factor(seurat_obj@active.ident, levels = pattern_order)

# Détection et filtrage des marqueurs
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
non_coding_pattern <- "^Gm|Rik|Mir|^mt-|^Rp[sl]|^Snr|^Malat1|^Xist"
markers_filtered <- markers %>%
  filter(!grepl(non_coding_pattern, gene, ignore.case = TRUE))

# Top 2 gènes par cluster
top3 <- markers_filtered %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

genes_top <- unique(top3$gene)
genes_top <- unique(c(genes_top[1:11], "Sftpa1", genes_top[12:length(genes_top)]))
genes_to_plot <- c(setdiff(genes_top, c("TdTomato", "Krt14")), "Krt14", "TdTomato")

# Définition des catégories par cluster
clusters_df <- data.frame(
  cluster = pattern_order,
  category = c(
    rep("Immune", 7),
    rep("Pulmonary", 2),
    rep("Cancerous", 3),
    rep("Unknown", 2)
  ),
  pos = seq_along(pattern_order)
)

category_df <- clusters_df %>%
  group_by(category) %>%
  summarise(
    ymin = min(pos) - 0.5,
    ymax = max(pos) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(color = case_when(
    category == "Immune" ~ "lightblue",
    category == "Pulmonary" ~ "lightgreen",
    category == "Cancerous" ~ "lightpink",
    category == "Unknown" ~ "lightgrey"
  ))

# Couleurs des clusters
pattern_colors <- c(
  "13" = "#ff66a8",
  "12" = "#fb61d7",
  "11" = "#df70f8",
  "10" = "#a58aff",
  "9"  = "#06a4ff",
  "8"  = "#00b6eb",
  "7"  = "#00bfc4",
  "6"  = "#00c094",
  "5"  = "#00bc56",
  "4"  = "#53b400",
  "3"  = "#99a800",
  "2"  = "#c49a00",
  "1"  = "#e38900",
  "0"  = "#f8766d"
)

colored_labels <- sapply(levels(seurat_obj@active.ident), function(x) {
  color <- pattern_colors[x]
  paste0("<span style='color:", color, "'>", x, "</span>")
})
names(colored_labels) <- levels(seurat_obj@active.ident)

labels_x <- sapply(genes_to_plot, function(gene) {
  if (gene %in% c("Krt14", "TdTomato")) {
    paste0("<span style='color:red; font-weight:bold'>", gene, "</span>")
  } else {
    gene
  }
})

# DotPlot
pA <- DotPlot(seurat_obj, features = genes_to_plot) +
  RotatedAxis() +
  scale_x_discrete(labels = labels_x) +
  scale_y_discrete(labels = colored_labels) +
  scale_color_gradient(low = "blue", high = "red") +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1, size = 6, face = "bold"),
    axis.text.y = element_markdown(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 50)
  ) +
  geom_rect(data = category_df,
            aes(xmin = -0.5, xmax = 0, ymin = ymin, ymax = ymax, fill = category),
            inherit.aes = FALSE, alpha = 0.3) +
  geom_text(data = category_df,
            aes(x = -0.25, y = (ymin + ymax) / 2, label = category),
            inherit.aes = FALSE,
            size = 3,
            fontface = "bold",
            angle = 90,
            color = ifelse(category_df$category == "Pulmonary", "darkgreen", "black")) +
  scale_fill_manual(
    values = c(
      "Immune" = "lightblue",
      "Pulmonary" = "lightgreen",
      "Cancerous" = "lightpink",
      "Unknown" = "lightgrey"
    ),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(-1, length(unique(genes_to_plot))))

# Chargement pour ImSig
seurat_obj2 <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds")
signature_imsig <- read.table("C:/Users/louis/Desktop/Stage/scrnaseq/signatures/ImSig2_mouse.txt", header = FALSE, stringsAsFactors = FALSE)$V1
seurat_obj2 <- AddModuleScore(seurat_obj2, features = list(signature_imsig), name = "ImSig_Score")

# Violin plot
pB <- VlnPlot(
  object = seurat_obj2,
  features = "ImSig_Score1",
  group.by = "seurat_clusters",
  pt.size = 0
) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
  labs(y = "ImSig Score", x = "Seurat clusters") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(color = "black", size = 16, margin = margin(t = 5)),
    axis.title.y = element_text(color = "black", size = 16),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

max_y1 <- max(seurat_obj2@meta.data$ImSig_Score1[!is.na(seurat_obj2@meta.data$seurat_clusters)], na.rm = TRUE)
pB <- pB + scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(NA, max_y1 * 1.15))

# Ajouter n
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
pB <- add_cell_count(pB, seurat_obj2@meta.data, "seurat_clusters")

# Panel combiné
panel_plot <- pA + pB + 
  plot_layout(ncol = 2, widths = c(1, 1)) + 
  plot_annotation(tag_levels = list(c("C", "D"))) & 
  theme(plot.tag = element_text(face = "bold", size = 20, hjust = 0, vjust = 1))

# Sauvegarde
ggsave("C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Rapport/Panel_ImSig_DotPlot_Seurat2.png",
       plot = panel_plot, width = 12, height = 6, dpi = 300)
