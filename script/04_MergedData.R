# Script: Merge and Quality Control Visualization of Multiple scRNA-seq Samples with Seurat
#
# Description:
# This script merges multiple filtered single-cell RNA sequencing (scRNA-seq) Seurat objects into a combined object.
# It performs dimensionality reduction (PCA, UMAP), clustering across multiple resolutions, and visualizes QC metrics.
# The script generates UMAP plots colored by various metadata such as counts, features, mitochondrial and ribosomal content,
# doublet scores, cell cycle phase, cluster identity, and sample origin.
# It also creates dot plots and cluster composition barplots by sample identity.
#
# Input:
# - Filtered Seurat RDS files located in `repertoire`
#
# Output:
# - QC and visualization plots saved as a single PDF file (`pdf_path`)
# - Combined Seurat object saved as an RDS file in `combine_sample_dir`
#
library(Seurat)
library(ggplot2)
library(clustree)
library(dplyr)

repertoire <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/03_FilteredData"
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/qc_plot_merged"
combine_sample_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData"
pdf_path <- file.path(output_dir, "QC_all_samples.pdf")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(combine_sample_dir, showWarnings = FALSE, recursive = TRUE)

fichiers_rds <- list.files(repertoire, pattern = "\\.rds$", full.names = TRUE)
objets_seurat <- lapply(fichiers_rds, readRDS) 

seurat_combined <- merge(objets_seurat[[1]], y = objets_seurat[-1], merge.data = TRUE)
seurat_combined[["RNA"]] <- JoinLayers(seurat_combined[["RNA"]])

cols_to_keep <- !grepl("^RNA_snn", colnames(seurat_combined@meta.data))
seurat_combined@meta.data <- seurat_combined@meta.data[, cols_to_keep]

seurat_combined <- RunPCA(seurat_combined, verbose = FALSE)

pdf(pdf_path, width = 10, height = 8)

seurat_combined <- FindNeighbors(seurat_combined)
seurat_combined <- FindClusters(seurat_combined, resolution = seq(0.1, 1, by = 0.1), verbose = FALSE)
clustree(seurat_combined, prefix = "RNA_snn_res.")

seurat_combined$seurat_clusters <- seurat_combined$RNA_snn_res.0.4
Idents(seurat_combined) <- seurat_combined$seurat_clusters

seurat_combined <- RunUMAP(seurat_combined, dims = 1:30)


num_dims <- length(seurat_combined@commands$RunUMAP$dims)
resolution <- "0.4" 


FeaturePlot(seurat_combined, features = "nCount_RNA", reduction = "umap") +
  scale_color_gradient(low = "black", high = "red") + 
  ggtitle(paste("UMAP Representation of nCount RNA\nDims:", num_dims, "- Resolution:", resolution))


FeaturePlot(seurat_combined, features = "nFeature_RNA", reduction = "umap") +
  scale_color_gradient(low = "black", high = "red") + 
  ggtitle(paste("UMAP Representation of nFeature RNA\nDims:", num_dims, "- Resolution:", resolution))

FeaturePlot(seurat_combined, features = "percent.mt", reduction = "umap") +
  scale_color_gradient(low = "black", high = "red") + 
  ggtitle(paste("UMAP Representation of Mitochondrial Percentage\nDims:", num_dims, "- Resolution:", resolution))

FeaturePlot(seurat_combined, features = "percent.ribo", reduction = "umap") +
  scale_color_gradient(low = "black", high = "red") + 
  ggtitle(paste("UMAP Representation of Ribosomal Percentage\nDims:", num_dims, "- Resolution:", resolution))

FeaturePlot(seurat_combined, features = "doublets.score", reduction = "umap") +
  scale_color_gradient(low = "black", high = "red") + 
  ggtitle(paste("UMAP Representation of Doublet Score\nDims:", num_dims, "- Resolution:", resolution))

DimPlot(seurat_combined, reduction = "umap", group.by = "Phase") +
  ggtitle(paste("UMAP Representation of Cell Cycle Phase\nDims:", num_dims, "- Resolution:", resolution))

DimPlot(seurat_combined, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle(paste("UMAP Representation of Clusters\nDims:", num_dims, "- Resolution:", resolution))

DimPlot(seurat_combined, reduction = "umap", group.by = "orig.ident") +
  ggtitle(paste("UMAP Representation by Sample Identity\nDims:", num_dims, "- Resolution:", resolution))

DotPlot(seurat_combined, features = "Vim", group.by = "orig.ident") +
  ggtitle("DotPlot of Vim Expression by Sample Identity")

FeaturePlot(seurat_combined, features = "Vim", reduction = "umap") +
  ggtitle(paste("UMAP Representation of Vim Expression\nDims:", num_dims, "- Resolution:", resolution))

FeaturePlot(seurat_combined, features = "TdTomato", reduction = "umap") +
  ggtitle(paste("UMAP Representation of TdTomato Expression\nDims:", num_dims, "- Resolution:", resolution))

FeaturePlot(seurat_combined, features = "SV40", reduction = "umap") +
  ggtitle(paste("UMAP Representation of SV40 Expression\nDims:", num_dims, "- Resolution:", resolution))

cluster_timepoint_df <- data.frame(
  Cluster = seurat_combined$seurat_clusters,
  TimePoint = seurat_combined$orig.ident
)


plot <- ggplot(cluster_timepoint_df, aes(x = Cluster, fill = TimePoint)) +
  geom_bar(position = "fill") +  
  theme_minimal(base_size = 15) +
  labs(
    title = "Time Point Distribution by Cluster (Proportions)",  
    x = "Cluster", 
    y = "Proportion"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  
    plot.title = element_text(hjust = 0.5),  
    panel.background = element_rect(fill = "white"),  
    plot.background = element_rect(fill = "white")   
  )

print(plot)

dev.off()

cat("✅ Fichier PDF enregistré :", pdf_path, "\n")

saveRDS(seurat_combined, file = file.path(combine_sample_dir, "seurat_combined.rds"))