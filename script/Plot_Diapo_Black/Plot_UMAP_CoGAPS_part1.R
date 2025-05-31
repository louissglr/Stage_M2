library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(forcats)
library(RColorBrewer)
library(cowplot)

# Lire l'objet Seurat combiné
seurat <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")

# Nettoyer les noms de clusters
seurat$pattern_cogaps <- gsub("Pattern_", "", seurat$pattern_cogaps)

# Ordre visuel 9 → 1 pour UMAP et barplot
cluster_order <- rev(sort(unique(seurat$pattern_cogaps)))
seurat$pattern_cogaps_visual <- factor(seurat$pattern_cogaps, levels = cluster_order)
seurat$pattern_cogaps_sizeorder <- factor(seurat$pattern_cogaps, levels = sort(unique(seurat$pattern_cogaps)))

# Préparer les données pour le barplot
df <- seurat@meta.data %>%
  mutate(pattern_cogaps_sizeorder = factor(pattern_cogaps_sizeorder, levels = levels(seurat$pattern_cogaps_sizeorder))) %>%
  count(pattern_cogaps_sizeorder, orig.ident, name = "n") %>%
  group_by(pattern_cogaps_sizeorder) %>%
  arrange(pattern_cogaps_sizeorder, desc(n)) %>%
  mutate(cluster_ident = paste0(pattern_cogaps_sizeorder, "_", orig.ident)) %>%
  ungroup()

df$cluster_ident <- factor(df$cluster_ident, levels = unique(df$cluster_ident))

# Ordre désiré pour la légende orig.ident
ordre_legende <- c("N1-0h", "N1-72h", "N1-1week", "N1-2week", "N1-4week", "N2-0h", "N2-72h", "N2-1week")

# Palette de couleurs Set2 assignée selon l’ordre
palette_orig <- setNames(brewer.pal(n = length(ordre_legende), "Set2"), ordre_legende)

# Couleurs pour chaque cluster_ident selon son orig.ident
cluster_ident_colors <- setNames(palette_orig[as.character(df$orig.ident)], df$cluster_ident)

# Palette UMAP couleurs pour clusters 9 → 1
umap_colors <- scales::hue_pal()(nlevels(seurat$pattern_cogaps_visual))
names(umap_colors) <- levels(seurat$pattern_cogaps_visual)

# Barplot (ordre des barres 9 → 1)
barplot <- ggplot(df, aes(x = pattern_cogaps_sizeorder, y = n, fill = cluster_ident)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_manual(values = cluster_ident_colors, guide = "none") +
  labs(x = "CoGAPS patterns", y = "Number of cells") +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text = element_text(size = 14, color = "white"),
    axis.text.y = element_text(color = umap_colors[as.character(levels(seurat$pattern_cogaps_sizeorder))], face = "bold"),
    axis.title = element_text(size = 15, color = "white"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white")
  )

# UMAP avec labels blancs
umap <- DimPlot(seurat, reduction = "umap", group.by = "pattern_cogaps_visual", label = TRUE, pt.size = 0.5, label.size = 8, label.color = "white") +
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


# Légende manuelle pour orig.ident
legend_data <- data.frame(orig.ident = ordre_legende)
legend_data$orig.ident <- factor(legend_data$orig.ident, levels = ordre_legende)

legend_plot <- ggplot(legend_data, aes(x = 1, y = orig.ident, color = orig.ident)) +
  geom_point(shape = 15, size = 8) +
  scale_color_manual(values = palette_orig, name = "Time points", breaks = ordre_legende) +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13, color = "white"),
    legend.text = element_text(size = 12, color = "white"),
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA)
  )

legend_only <- cowplot::get_legend(legend_plot)
legend_gg <- cowplot::ggdraw(legend_only)

# Assemblage final sans tags
panel <- umap + barplot + legend_gg +
  plot_layout(widths = c(1.3, 1.1, 0.3)) +
  plot_annotation() &  # PAS de tag (vide)
  theme(
    plot.background = element_rect(fill = "black", color = NA)
  )

# Sauvegarde
ggsave(
  filename = "C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Diapo_Black/cogaps_UMAP.png",
  plot = panel,
  width = 13,
  height = 6,
  dpi = 300,
  bg = "black"
)
