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
