library(Seurat)
library(future)
library(ggplot2)

# Augmenter la mémoire globale pour les objets futurs
options(future.globals.maxSize = 30 * 1024^3)  # 30 Go

# 1. Lire tous les fichiers RDS
files <- list.files("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/03_FilteredData", pattern = "\\.rds$", full.names = TRUE)

seurat_list <- lapply(files, readRDS)

# 1bis. Ajouter l'objet du Pérou
perou_rds <- readRDS("C:/Users/louis/Desktop/Stage/perou/output/04-IntegratedData/seurat_integrated_filtered.rds")
perou_rds$orig.ident <- "Ref_Perou"
seurat_list <- c(seurat_list, list(perou_rds))

# 2. Fusionner les objets Seurat
seurat_merged <- Reduce(function(x, y) merge(x, y), seurat_list)
seurat_merged$orig.ident <- gsub("SeuratProject", "Ref_Perou", seurat_merged$orig.ident)
seurat_merged <- SetIdent(seurat_merged, value = "orig.ident")
seurat_merged <- JoinLayers(seurat_merged)

seurat_merged <- NormalizeData(seurat_merged)
hvgs_old <- VariableFeatures(seurat_merged)

# 7. Split final par identifiant d’origine
seurat_merged[["RNA"]] <- split(seurat_merged[["RNA"]], f = seurat_merged$orig.ident)

# Détection des HVGs
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Récupérer les HVGs pour chaque couche
data.layers <- Layers(seurat_merged)[grep("data.", Layers(seurat_merged))]
print(data.layers)

hvgs_per_dataset <- lapply(data.layers, function(x) VariableFeatures(seurat_merged, layer = x))
names(hvgs_per_dataset) <- data.layers

# Ajouter les HVGs globaux et anciens
hvgs_per_dataset$all <- VariableFeatures(seurat_merged)
hvgs_per_dataset$old <- hvgs_old

# Calcul de l’overlap
temp <- unique(unlist(hvgs_per_dataset))
overlap <- sapply(hvgs_per_dataset, function(x) { temp %in% x })
pheatmap::pheatmap(t(overlap * 1), cluster_rows = FALSE,
                   color = c("grey90", "grey20"))
hvgs_all <- hvgs_per_dataset$all

# Scaling, PCA
seurat_merged <- ScaleData(seurat_merged, features = hvgs_all)
seurat_merged <- RunPCA(seurat_merged, features = hvgs_all, verbose = FALSE)

# Intégration CCA
seurat_merged <- IntegrateLayers(
  object = seurat_merged,
  method = CCAIntegration, orig.reduction = "pca",
  new.reduction = "integrated_cca", verbose = FALSE
)
seurat_merged <- FindNeighbors(seurat_merged, reduction = "integrated_cca", dims = 1:30)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.4, cluster.name = "cca_clusters")
seurat_merged <- RunUMAP(seurat_merged, reduction = "integrated_cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(seurat_merged, reduction = "umap.cca", label = TRUE)

# Intégration RPCA
seurat_merged <- IntegrateLayers(
  object = seurat_merged,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated_rpca",
  k.weight = 50,   # <- Abaisse à 50 (ou 30, 20 selon ton dataset)
  verbose = FALSE
)

seurat_merged <- FindNeighbors(seurat_merged, reduction = "integrated_rpca", dims = 1:30)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.4, cluster.name = "rpca_clusters")
seurat_merged <- RunUMAP(seurat_merged, reduction = "integrated_rpca", dims = 1:30, reduction.name = "umap.rpca")
DimPlot(seurat_merged, reduction = "umap.rpca", label = TRUE)

# Intégration Harmony
seurat_merged <- IntegrateLayers(
  object = seurat_merged, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
seurat_merged <- FindNeighbors(seurat_merged, reduction = "harmony", dims = 1:30)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.4, cluster.name = "harmony_clusters")
seurat_merged <- RunUMAP(seurat_merged, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(seurat_merged, reduction = "umap.harmony", label = TRUE)

# Vérification des réductions
names(seurat_merged@reductions)

# Sauvegarde finale
saveRDS(seurat_merged, file = "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/15_IntegrationWithPerou/integration.rds")
