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
  
  cols <- c("Sig" = "red", "Not Sig" = "grey")
  
  p <- ggplot(cluster_data, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300))) +
    geom_point(aes(color = signif), alpha = 0.9, size = 1.8, show.legend = show_legend) +
    scale_color_manual(values = cols, name = "Significance") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(
      title = paste("Cluster", clust_id),
      x = expression(Log[2]~"fold change"),
      y = expression(-Log[10]~"adjusted p-value")
    ) +
    xlim(xlim_vals) +
    ylim(ylim_vals) +
    
    # THÈME CLAIR
    theme_minimal(base_family = "sans") +
    theme(
      legend.position = if (show_legend) "right" else "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold")
    ) +
    
    # Labels top gènes
    geom_text_repel(
      data = subset(cluster_data, !is.na(label)),
      aes(label = label),
      color = "black",
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
      color = "black",
      size = 4,
      hjust = 0
    )
  }
  
  return(p)
}

# Génération des trois volcano plots
p2 <- volcano_plot_fixed(markers_all, 2, x_range, y_range, show_legend = FALSE)
p3 <- volcano_plot_fixed(markers_all, 3, x_range, y_range, show_legend = FALSE)
p9 <- volcano_plot_fixed(markers_all, 9, x_range, y_range, show_legend = TRUE)

# Combinaison et export
combined_plot <- p2 + p3 + p9 + plot_layout(ncol = 3)
print(combined_plot)

ggsave("volcano_clusters_2_3_9_top10_axes.png", combined_plot, width = 12, height = 4, dpi = 300)
