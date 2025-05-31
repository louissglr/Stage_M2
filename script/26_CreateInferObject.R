# Chargement de l'objet Seurat
seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/21_AddCoGAPSToSeuratCancer/seurat_cogaps_res_6p.rds")
seurat <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds")
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


# Étape 1 : Sélection des cellules "immune" selon CoGAPS_res_9p
immune_cells <- colnames(seurat)[seurat$CoGAPS_res_9p %in% c("Pattern_1", "Pattern_4", "Pattern_5")]

# Étape 2 : Sous-ensemble de ces cellules
seurat_immune <- subset(seurat, cells = immune_cells)

# Étape 3 : Annotation simple en "immune"
seurat_immune$pattern_cogaps_immune <- "immune"

# Étape 4 : Copie des annotations dans seurat_obj
seurat_obj$pattern_cogaps_immune <- seurat_obj$pattern_cogaps

# Étape 5 : Fusion propre
merged_obj <- merge(seurat_obj, y = seurat_immune, add.cell.ids = c("Main", "Immune"), merge.data = TRUE)

# Étape 6 : Vérification
table(merged_obj$pattern_cogaps_immune)




# Noms des cellules
cells_all <- colnames(merged_obj)

# Identifier les cellules immune (pattern_cogaps_immune == "immune")
immune_cells <- cells_all[merged_obj$pattern_cogaps_immune == "immune"]

# Identifier cellules N1 non-immune
cells_N1_nonimmune <- grep("^Main_N1", cells_all, value = TRUE)
cells_N1_nonimmune <- setdiff(cells_N1_nonimmune, immune_cells)

# Identifier cellules N2 non-immune
cells_N2_nonimmune <- grep("^Main_N2", cells_all, value = TRUE)
cells_N2_nonimmune <- setdiff(cells_N2_nonimmune, immune_cells)

# Créer le vecteur cellules pour chaque objet en ajoutant TOUTES les immune
cells_N1 <- c(cells_N1_nonimmune, immune_cells)
cells_N2 <- c(cells_N2_nonimmune, immune_cells)

# Créer les objets subset
merged_obj_N1 <- subset(merged_obj, cells = cells_N1)
merged_obj_N2 <- subset(merged_obj, cells = cells_N2)

# Vérification rapide
cat("Table pattern_cogaps_immune - merged_obj_N1\n")
print(table(merged_obj_N1$pattern_cogaps_immune))

cat("Table pattern_cogaps_immune - merged_obj_N2\n")
print(table(merged_obj_N2$pattern_cogaps_immune))

# Sauvegarde
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/26_CreateInferObjectCancer"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveRDS(merged_obj_N1, file = file.path(output_dir, "merged_obj_N1.rds"))
saveRDS(merged_obj_N2, file = file.path(output_dir, "merged_obj_N2.rds"))
saveRDS(merged_obj, file = file.path(output_dir, "merged_obj.rds"))
