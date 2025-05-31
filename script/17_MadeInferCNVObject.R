# Chargement de l'objet Seurat
seurat <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds")

# Étape 1 : annotation avec les clusters Seurat
seurat$annot <- seurat$seurat_clusters

# Étape 2 : annotation des patterns CoGAPS
res_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/08_CoGAPS"
rds_files <- list.files(res_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

for (file_path in rds_files) {
  cogaps <- readRDS(file_path)
  patterns_in_order <- t(cogaps@sampleFactors[colnames(seurat), , drop = FALSE])
  dominant_patterns <- apply(patterns_in_order, 2, function(x) {
    rownames(patterns_in_order)[which.max(x)]
  })
  col_name <- paste0("CoGAPS_", basename(dirname(file_path)))
  seurat[[col_name]] <- dominant_patterns
}

# Récupération des valeurs pour filtrage
patterns <- seurat$CoGAPS_res_9p
idents <- seurat$orig.ident

# Ajout de la colonne group_label à l'objet total
seurat$group_label <- ifelse(
  seurat$CoGAPS_res_9p %in% c("Pattern_4", "Pattern_5"),
  seurat$CoGAPS_res_9p,
  seurat$orig.ident
)

# Création de seurat_obj1 : N1 ou N2 + Pattern_4/5
cells_obj1 <- colnames(seurat)[
  grepl("^N1", idents) |
    (grepl("^N2", idents) & patterns %in% c("Pattern_4", "Pattern_5"))
]
seurat_obj1 <- subset(seurat, cells = cells_obj1)

# Ajout de la colonne group_label à seurat_obj1
seurat_obj1$group_label <- ifelse(
  seurat_obj1$CoGAPS_res_9p %in% c("Pattern_4", "Pattern_5"),
  seurat_obj1$CoGAPS_res_9p,
  seurat_obj1$orig.ident
)

# Création de seurat_obj2 : N2 ou N1 + Pattern_4/5
cells_obj2 <- colnames(seurat)[
  grepl("^N2", idents) |
    (grepl("^N1", idents) & patterns %in% c("Pattern_4", "Pattern_5"))
]
seurat_obj2 <- subset(seurat, cells = cells_obj2)

# Ajout de la colonne group_label à seurat_obj2
seurat_obj2$group_label <- ifelse(
  seurat_obj2$CoGAPS_res_9p %in% c("Pattern_4", "Pattern_5"),
  seurat_obj2$CoGAPS_res_9p,
  seurat_obj2$orig.ident
)

# Création du dossier de sortie si nécessaire
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/17_MadeInferCNVObject"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Sauvegarde des objets
saveRDS(seurat_obj1, file = file.path(output_dir, "seurat_obj1.rds"))
saveRDS(seurat_obj2, file = file.path(output_dir, "seurat_obj2.rds"))
saveRDS(seurat, file = file.path(output_dir, "seurat_cogaps_res_9p.rds"))