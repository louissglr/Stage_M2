library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)

# Load your Seurat object
seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds")

# Step 1: Extract metadata and calculate the proportion for each cluster and each orig.ident
df <- seurat_obj@meta.data %>%
  count(seurat_clusters, orig.ident) %>%
  group_by(seurat_clusters) %>%
  mutate(
    prop = n / sum(n)  # Calculate the proportion of each orig.ident in each cluster
  ) %>%
  ungroup()

# Step 2: Barplot with black background, but without labels inside the bars
barplot <- ggplot(df, aes(x = factor(seurat_clusters), y = n, fill = factor(orig.ident))) +
  geom_bar(stat = "identity", position = "stack", color = "white") +  # Stacked bars
  labs(
    title = "Cell Proportion by Cluster and Timepoint (orig.ident)",
    x = "Cluster",
    y = "Number of Cells",
    fill = "Timepoint (orig.ident)"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(color = "white", size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "white", size = 16),
    axis.title.x = element_text(color = "white", size = 18, face = "bold"),
    axis.title.y = element_text(color = "white", size = 18, face = "bold"),
    legend.text = element_text(color = "white", size = 14),
    legend.title = element_text(color = "white", size = 14, face = "bold"),
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"),
    plot.title = element_text(hjust = 0.5, color = "white", size = 18, face = "bold")  # Increased title size
  )

# Display the plot
print(barplot)

# Save the plot as a PNG
ggsave(
  "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation/barplot_timepoint_clusters_no_labels.png", 
  plot = barplot, 
  width = 10,    # Image width in inches
  height = 8,    # Image height in inches
  dpi = 300      # Image resolution
)

# Créer le dataframe avec les colonnes souhaitées
df <- seurat_obj@meta.data %>%
  count(seurat_clusters, orig.ident) %>%
  rename(
    clusters = seurat_clusters,
    timepoint = orig.ident,
    nb_cellules = n
  ) %>%
  arrange(clusters, timepoint)

# Vérification rapide des premières lignes du dataframe
print(head(df))

# Sauvegarder le dataframe dans un fichier CSV
write_csv(df, "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/UMAP_presentation/nb_cellules_par_cluster_et_timepoint.csv")
