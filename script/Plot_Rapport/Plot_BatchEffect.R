library(ggplot2)
library(dplyr)
library(scales)
library(openxlsx)
library(VennDiagram)
library(grid)
library(patchwork)
library(RColorBrewer)

# --- Chargement des données cogaps ---
cogaps <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")
cogaps$pattern_cogaps <- gsub("Pattern_", "", cogaps$pattern_cogaps)

# --- Clusters à inclure artificiellement par pattern ---
clusters_to_add <- tibble::tribble(
  ~pattern_cogaps, ~seurat_clusters,
  "4", 6,
  "4", 11,
  "4", 12,
  "4", 13,
  "5", 6,
  "5", 9,
  "5", 11,
  "5", 13
)

# --- Dataframe pour barplot 1 (pattern 4 et 5) ---
df_pattern_init <- cogaps@meta.data %>%
  filter(pattern_cogaps %in% c("4", "5")) %>%
  count(pattern_cogaps, seurat_clusters) %>%
  mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))

missing_rows <- clusters_to_add %>%
  anti_join(df_pattern_init, by = c("pattern_cogaps", "seurat_clusters")) %>%
  mutate(n = 0)

df_pattern <- bind_rows(df_pattern_init, missing_rows) %>%
  group_by(pattern_cogaps) %>%
  mutate(
    prop = n / sum(n),
    label = ifelse(prop >= 0.15, as.character(seurat_clusters), NA)
  ) %>%
  # ORDRE PAR TAILLE DÉCROISSANTE dans chaque pattern pour empilement correct
  arrange(pattern_cogaps, desc(n)) %>%
  mutate(
    pos_n = cumsum(n) - n / 2,
    rang = row_number()  # rang dans chaque pattern selon taille décroissante
  ) %>%
  ungroup()

df_pattern$seurat_clusters <- factor(df_pattern$seurat_clusters)
totals_pattern <- df_pattern %>%
  group_by(pattern_cogaps) %>%
  summarise(total_cells = sum(n))

cluster_levels <- levels(df_pattern$seurat_clusters)
palette_colors <- scales::hue_pal()(length(cluster_levels))
names(palette_colors) <- cluster_levels

barplot1 <- ggplot(df_pattern, aes(x = factor(pattern_cogaps), y = n, fill = seurat_clusters,
                                   group = interaction(pattern_cogaps, -rang))) +
  geom_bar(stat = "identity", color = "white") +
  geom_text(aes(y = pos_n, label = label), color = "black", size = 5, fontface = "bold", na.rm = TRUE) +
  geom_text(data = totals_pattern, aes(x = pattern_cogaps, y = total_cells + 50,
                                       label = paste0("n = ", total_cells)),
            inherit.aes = FALSE, size = 3, hjust = 0, fontface = "plain", color = "black") +
  coord_flip(clip = "off") +
  ylim(0, 2300) +
  labs(x = "Pattern", y = "Number of Cells", fill = "Cluster") +
  scale_fill_manual(values = palette_colors) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(color = "black", size = 16),
    axis.text.y = element_text(color = "black", size = 18, face = "bold"),
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 14),
    legend.title = element_text(color = "black", size = 14),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# --- Barplot 2 ---

df_cluster <- cogaps@meta.data %>%
  filter(seurat_clusters %in% c("2", "3", "9")) %>%
  count(seurat_clusters, orig.ident) %>%
  group_by(seurat_clusters) %>%
  mutate(
    prop = n / sum(n),
    label = ifelse(prop >= 0.2, as.character(orig.ident), NA)
  ) %>%
  arrange(seurat_clusters, desc(n)) %>%
  mutate(
    pos_n = cumsum(n) - n / 2,
    rang = row_number()
  ) %>%
  ungroup()

totals_cluster <- df_cluster %>%
  group_by(seurat_clusters) %>%
  summarise(total_cells = sum(n))

ordre_legende <- c("N1-0h", "N1-72h", "N1-1week", "N1-2week", "N1-4week", "N2-0h", "N2-72h", "N2-1week")
df_cluster$orig.ident <- factor(df_cluster$orig.ident, levels = ordre_legende)
palette_orig <- setNames(brewer.pal(n = length(ordre_legende), "Set2"), ordre_legende)

cluster_colors <- c("2" = "#c49a00", "3" = "#99a800", "9" = "#06a4ff")
df_cluster$seurat_clusters <- factor(df_cluster$seurat_clusters, levels = names(cluster_colors))

barplot2 <- ggplot(df_cluster, aes(x = factor(seurat_clusters), y = n, fill = orig.ident,
                                   group = interaction(seurat_clusters, -rang))) +
  geom_bar(stat = "identity", color = "white") +
  geom_text(aes(y = pos_n, label = label), color = "black", size = 5, fontface = "bold", na.rm = TRUE) +
  geom_text(data = totals_cluster, aes(x = seurat_clusters, y = total_cells + 50,
                                       label = paste0("n = ", total_cells)),
            inherit.aes = FALSE, size = 3, hjust = 0, fontface = "plain", color = "black") +
  coord_flip(clip = "off") +
  labs(x = "Seurat Cluster", y = "Number of Cells", fill = "Time points") +
  scale_fill_manual(values = palette_orig) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(color = "black", size = 16),
    axis.text.y = element_text(
      size = 18, face = "bold",
      color = cluster_colors[levels(df_cluster$seurat_clusters)]
    ),
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 14),
    legend.title = element_text(color = "black", size = 14),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# --- Venn Diagram ---
file_path <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/06_FindAllMarkers/Markers_seurat_combined.rds.xlsx"
markers_all <- read.xlsx(file_path, sheet = "res_0.4_significant")

clusters_of_interest <- c(2, 3, 9)
filtered <- markers_all %>% filter(cluster %in% clusters_of_interest)
filtered_genes <- filtered %>%
  filter(pct.1 >= 0.25, avg_log2FC >= 0.25)

gene_lists <- split(filtered_genes$gene, filtered_genes$cluster)
names(gene_lists) <- paste0("Cluster ", names(gene_lists))

cluster_order <- c("Cluster 2", "Cluster 3", "Cluster 9")
gene_lists <- gene_lists[cluster_order]
colors <- c("#c49a00", "#99a800", "#06a4ff")

venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = cluster_order,
  filename = NULL,
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 500,
  fill = colors,
  col = "black",
  cex = 2.2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 2.0,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.col = colors,
  cat.dist = 0.06,
  cat.pos = c(-10, 10, -150),
  margin = 0.05,
  bg = "white"
)

venn_grob <- grobTree(venn.plot)
venn_gg <- ggplot() + annotation_custom(venn_grob) + theme_void()

# --- Patchwork final ---
combined_plot <- barplot1 + venn_gg + barplot2 + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(size = 24, face = "bold", family = "Arial"),
    plot.tag.position = c(0.02, 0.98),
    legend.position = "right"
  )

# --- Export PNG ---
ggsave("C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Rapport/panel_ABC_patchwork.png",
       combined_plot, width = 18, height = 6, dpi = 300)
