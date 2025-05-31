library(openxlsx)
library(VennDiagram)
library(dplyr)

# Charger le fichier Excel
file_path <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/06_FindAllMarkers/Markers_seurat_combined.rds.xlsx"
markers_all <- read.xlsx(file_path, sheet = "res_0.4_significant")

# Filtrer pour les clusters d'intérêt
clusters_of_interest <- c(2, 3, 9)
filtered <- markers_all %>% filter(cluster %in% clusters_of_interest)

# Apply filtering based on min.pct and logfc.threshold
filtered_genes <- filtered %>%
  filter(pct.1 >= 0.25, avg_log2FC >= 0.25)

# Extraire les gènes pour chaque cluster
gene_lists <- split(filtered_genes$gene, filtered_genes$cluster)

# Renommer les listes pour le Venn
names(gene_lists) <- paste0("Cluster_", names(gene_lists))

# Définir le chemin de sortie
output_path <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation/venn_clusters.png"

# Ouvrir le fichier PNG avec fond noir
png(output_path, width = 3000, height = 3000, res = 500, bg = "black")

# Créer le diagramme de Venn sans l’enregistrer tout de suite
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = names(gene_lists),
  filename = NULL,
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 500,
  fill = c("red", "blue", "green"),
  col = "white",                 # Contour blanc
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.8,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.col = "white",            # Étiquettes blanches
  cat.dist = 0.07,
  cat.pos = c(-10, 10, -150),    # Position des étiquettes (ex : Cluster9 en bas)
  margin = 0.1
)

# Dessiner le fond noir
grid.newpage()
grid.rect(gp = gpar(fill = "black", col = NA))
grid.draw(venn.plot)

# Sauvegarder
dev.off()



library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(openxlsx)
file_path <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/06_FindAllMarkers/Markers_seurat_combined.rds.xlsx"
markers_all <- read.xlsx(file_path, sheet = "res_0.4_all")


# 1. Calcul des bornes globales
markers_subset <- markers_all %>% filter(cluster %in% c(2, 3, 9))
markers_subset <- markers_subset %>%
  filter(pct.1 >= 0.25, avg_log2FC >= 0.25)

x_range <- range(markers_subset$avg_log2FC, na.rm = TRUE)
y_range <- range(-log10(markers_subset$p_val_adj + 1e-310), na.rm = TRUE)  # éviter log10(0)

volcano_plot_fixed <- function(data, clust_id, xlim_vals, ylim_vals, show_legend = FALSE) {
  cluster_data <- data %>% filter(cluster == clust_id)
  
  cluster_data <- cluster_data %>%
    mutate(
      signif = ifelse(p_val_adj < 0.05, "Sig", "Not Sig")
    )
  
  top10 <- cluster_data %>%
    filter(signif == "Sig") %>%
    arrange(p_val_adj) %>%
    slice_head(n = 10)
  
  cluster_data$label <- ifelse(cluster_data$gene %in% top10$gene, cluster_data$gene, NA)
  
  cols <- c("Sig" = "red", "Not Sig" = "white")
  
  p <- ggplot(cluster_data, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300))) +
    geom_point(aes(color = signif), alpha = 0.9, size = 1.8, show.legend = show_legend) +
    scale_color_manual(values = cols, name = "Significance") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "white") +
    labs(
      title = paste("Cluster", clust_id),
      x = expression(Log[2]~"fold change"),
      y = expression(-Log[10]~"adjusted p-value")
    ) +
    xlim(xlim_vals) +
    ylim(ylim_vals) +
    
    # THÈME NOIR
    theme_minimal(base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      legend.background = element_rect(fill = "black", color = NA),
      legend.box.background = element_rect(fill = "black", color = NA),
      panel.grid = element_blank(),
      
      axis.text = element_text(color = "white", size = 12),
      axis.title = element_text(color = "white", size = 14, face = "bold"),
      axis.line = element_line(color = "white"),
      axis.ticks = element_line(color = "white"),
      
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(color = "white", size = 12),
      legend.title = element_text(color = "white", size = 14, face = "bold"),
      
      plot.title = element_text(hjust = 0.5, color = "white", size = 16, face = "bold")
    ) +
    
    # Labels top gènes
    geom_text_repel(
      data = subset(cluster_data, !is.na(label)),
      aes(label = label),
      color = "white",
      size = 3,
      fontface = "italic",
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "gray70",
      segment.size = 0.3,
      min.segment.length = 0.1,
      force = 2,
      max.overlaps = Inf,
      seed = 42
    )
  
  # Ajouter le label p=0.05 uniquement pour le cluster 9
  if (clust_id == 9) {
    p <- p + annotate(
      "text",
      x = max(xlim_vals) * 0.95,
      y = -log10(0.05),
      label = "p = 0.05",
      color = "white",
      size = 4,
      hjust = 0
    )
  }
  
  return(p)
}




p2 <- volcano_plot_fixed(markers_all, 2, x_range, y_range, show_legend = FALSE)
p3 <- volcano_plot_fixed(markers_all, 3, x_range, y_range, show_legend = FALSE)
p9 <- volcano_plot_fixed(markers_all, 9, x_range, y_range, show_legend = TRUE)

combined_plot <- p2 + p3 + p9 + plot_layout(ncol = 3)
print(combined_plot)

ggsave("volcano_clusters_2_3_9_top10_axes.png", combined_plot, width = 12, height = 4, dpi = 300)
