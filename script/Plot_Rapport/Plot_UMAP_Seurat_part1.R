library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(forcats)
library(RColorBrewer)
library(cowplot)

# Charger l'objet Seurat
seurat <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds")

# Ordonner les clusters par taille décroissante (pour barplot)
cluster_sizes <- seurat@meta.data %>%
  count(seurat_clusters, name = "total_cells") %>%
  arrange(desc(total_cells))

# Fixer l'ordre des clusters pour le barplot (par taille décroissante)
seurat$seurat_clusters <- factor(seurat$seurat_clusters, levels = cluster_sizes$seurat_clusters)

# Préparer les données pour le barplot
df <- seurat@meta.data %>%
  count(seurat_clusters, orig.ident, name = "n") %>%
  group_by(seurat_clusters) %>%
  arrange(seurat_clusters, desc(n)) %>%
  mutate(cluster_ident = paste0(seurat_clusters, "_", orig.ident)) %>%
  ungroup()

df$cluster_ident <- factor(df$cluster_ident, levels = unique(df$cluster_ident))

# Ordre désiré pour la légende orig.ident
ordre_legende <- c("N1-0h", "N1-72h", "N1-1week", "N1-2week", "N1-4week", "N2-0h", "N2-72h", "N2-1week")

# Palette de couleurs Set2 assignée selon l’ordre
palette_orig <- setNames(brewer.pal(n = length(ordre_legende), "Set2"), ordre_legende)

# Couleurs pour chaque cluster_ident selon son orig.ident
cluster_ident_colors <- setNames(palette_orig[as.character(df$orig.ident)], df$cluster_ident)

# Palette UMAP couleurs pour clusters, dans l’ordre naturel (0,1,2,...)
umap_colors <- scales::hue_pal()(nlevels(seurat$seurat_clusters))
names(umap_colors) <- levels(seurat$seurat_clusters)

# Barplot avec labels en gras et colorés selon la palette umap_colors (en respectant l’ordre des clusters du barplot)
barplot <- ggplot(df, aes(x = seurat_clusters, y = n, fill = cluster_ident)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_manual(values = cluster_ident_colors, guide = "none") +
  labs(x = "Seurat clusters", y = "Number of cells") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.text.y = element_text(color = umap_colors[levels(seurat$seurat_clusters)], face = "bold"),
    axis.title = element_text(size = 15, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# UMAP avec couleurs correctes des clusters (ordre naturel)
umap <- DimPlot(seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
  scale_color_manual(values = umap_colors, name = "Clusters") +
  # Ici on inverse juste l'ordre de la légende avec guide_legend, mais sans modifier les données
  guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 6))) +
  labs(title = NULL) +
  theme(
    plot.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# Création de la légende manuelle pour orig.ident
legend_data <- data.frame(orig.ident = ordre_legende)
legend_data$orig.ident <- factor(legend_data$orig.ident, levels = ordre_legende)

legend_plot <- ggplot(legend_data, aes(x = 1, y = orig.ident, color = orig.ident)) +
  geom_point(shape = 15, size = 8) +
  scale_color_manual(values = palette_orig, name = "Time points", breaks = ordre_legende) +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.background = element_rect(fill = "white", color = NA)
  )

legend_only <- cowplot::get_legend(legend_plot)
legend_gg <- cowplot::ggdraw(legend_only)

# Assemblage final
panel <- umap + barplot + legend_gg +
  plot_layout(widths = c(1.3, 1.1, 0.3)) +
  plot_annotation(tag_levels = list(c("A", "B", ""))) &  
  theme(
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Sauvegarde
ggsave(
  filename = "C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Rapport/seurat_UMAP.png",
  plot = panel,
  width = 13,
  height = 6,
  dpi = 300
)
