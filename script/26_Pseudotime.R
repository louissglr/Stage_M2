library(Seurat)
library(monocle3)
library(ggplot2)
library(patchwork)
library(dplyr)
library(forcats)
library(RColorBrewer)
library(cowplot)
library(SeuratWrappers)

# Charger l'objet Seurat
seurat <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/21_AddCoGAPSToSeuratCancer/seurat_cogaps_res_6p.rds")
DefaultAssay(seurat) <- "RNA"
# Nettoyer les noms de patterns
seurat$pattern_cogaps <- gsub("Pattern_", "", seurat$pattern_cogaps)

# Ordre visuel 9 → 1 pour UMAP et barplot
cluster_order <- rev(sort(unique(seurat$pattern_cogaps)))
seurat$pattern_cogaps_visual <- factor(seurat$pattern_cogaps, levels = cluster_order)

# Palette UMAP couleurs pour clusters 9 → 1 (teinte par défaut hue)
umap_colors <- scales::hue_pal()(nlevels(seurat$pattern_cogaps_visual))
names(umap_colors) <- levels(seurat$pattern_cogaps_visual)

# UMAP patterns (Seurat) avec légende organisée, labels blancs, couleurs personnalisées
# UMAP patterns (Seurat) avec légende organisée, labels blancs, couleurs personnalisées
p_left <- DimPlot(seurat, reduction = "umap", group.by = "pattern_cogaps_visual", label = TRUE, pt.size = 1.1, label.size = 8, label.color = "white") +
  scale_color_manual(values = umap_colors, name = "Patterns") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(title = NULL) +
  theme(
    plot.title = element_blank(),
    legend.text = element_text(size = 12, color = "white"),
    legend.title = element_text(size = 14, color = "white"),
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


# Convertir en CellDataSet Monocle3
cds <- as.cell_data_set(seurat)
# Ajouter les noms de gènes dans rowData pour l'affichage
rowData(cds)$gene_short_name <- rownames(cds)

cds@colData$pattern_cogaps <- seurat$pattern_cogaps
cds@int_colData$reducedDims$UMAP <- seurat[["umap"]]@cell.embeddings

# Trajectoire Monocle3
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
root_cell <- colnames(seurat)[seurat$pattern_cogaps == "3"][8]
cds <- order_cells(cds, root_cells = root_cell)
root_coords <- reducedDims(cds)$UMAP[root_cell, ]

# Plot UMAP pseudotime + trajectoire (Monocle3) avec taille des cellules augmentée
p_right <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  trajectory_graph_color = "white",
  cell_size = 1.2  # Taille augmentée des cellules
) +
  geom_point(aes(x = root_coords[1], y = root_coords[2]),
             color = "red", size = 3, shape = 4, stroke = 2) +
  geom_text(aes(x = root_coords[1], y = root_coords[2], label = "Start"),
            vjust = -1, color = "red", size = 4) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid = element_blank(),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white", size = 14),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"),
    plot.title = element_text(color = "white", size = 16, hjust = 0.5),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  ) +
  ggtitle("UMAP - Pseudotime")

# Combiner plots côte à côte avec patchwork
combined_plot <- p_left + p_right + plot_layout(ncol = 2)

# Afficher
print(combined_plot)

# Sauvegarder
ggsave(
  filename = "C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Diapo_Black/pseudotime_6p.png",
  plot = combined_plot,
  width = 16,
  height = 7,
  dpi = 300,
  bg = "black"
)

deg_pseudotime <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# Sélectionner les gènes significativement associés au pseudotime
deg_sig <- deg_pseudotime %>%
  filter(q_value < 0.05) %>%
  arrange(q_value)


# Filtrer les gènes significatifs si ce n’est pas déjà fait
deg_sig <- deg_pseudotime %>%
  filter(q_value < 0.05) %>%
  arrange(q_value)

# Sauvegarder dans un CSV
write.csv(deg_sig,
          file = "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/25_Pseudotime/deg_pseudotime_significant.csv",
          row.names = FALSE)



# Assure-toi que les noms courts sont bien renseignés
rowData(cds)$gene_short_name <- rownames(cds)

# Liste des gènes à tracer
genes_to_plot <- c("Vim", "Krt14", "TdTomato")

# Violin plot par pattern CoGAPS
p_violin <- plot_genes_violin(
  cds[genes_to_plot, ],
  group_cells_by = "pattern_cogaps",
  label_by_short_name = TRUE
)

# Optionnel : ajouter ta palette personnalisée si tu veux unifier les couleurs avec UMAP
p_violin <- p_violin + scale_fill_manual(values = umap_colors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Afficher
print(p_violin)











# 9. Visualiser quelques gènes significatifs le long du pseudotemps
plot_genes_in_pseudotime(cds[row.names(deg_sig)[1:20], ], color_cells_by = "pseudotime")


# Identifier le rowname correspondant à Vim
vim_row <- rownames(cds)[rowData(cds)$gene_short_name == "Krt14"]

# Visualiser Vim en fonction du pseudotemps
plot_genes_in_pseudotime(cds[vim_row, ], color_cells_by = "pseudotime")
