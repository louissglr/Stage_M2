library(Seurat)
library(CoGAPS)

# Charger l'objet Seurat avec résultats CoGAPS
cogaps <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")

# S'assurer que l'assay RNA est actif
DefaultAssay(cogaps) <- "RNA"

# Avant filtrage
nb_genes_avant <- nrow(cogaps[["RNA"]])
cat("Nombre de gènes avant filtrage :", nb_genes_avant, "\n")

# Filtrer les gènes qui ne commencent PAS par Gm, Rik, ou LOC
genes_to_keep <- !grepl("^Gm|^Rik|^LOC", rownames(cogaps[["RNA"]]), ignore.case = FALSE)

# Afficher le nombre de gènes concernés par le filtrage
nb_genes_retires <- sum(!genes_to_keep)
nb_genes_restants <- sum(genes_to_keep)

cat("Nombre de gènes retirés :", nb_genes_retires, "\n")
cat("Nombre de gènes conservés :", nb_genes_restants, "\n")

# Appliquer le filtrage
cogaps <- cogaps[genes_to_keep, ]

# Subsetting des cellules avec patterns 2, 3, 7, 8, 9
patterns_to_keep <- paste0("Pattern_", c(2, 3, 7, 8, 9))
cogaps_subset <- subset(cogaps, subset = pattern_cogaps %in% patterns_to_keep)

# Vérifier le résultat
cat("Distribution des cellules par pattern :\n")
print(table(cogaps_subset$pattern_cogaps))

cat("Distribution des cellules par origine :\n")
print(table(cogaps_subset$orig.ident))

# Définir le chemin du dossier de sortie
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/19_Keep_CoGAPS_Cancer"

# Créer le dossier s'il n'existe pas
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Définir le chemin complet pour le fichier RDS
output_file <- file.path(output_dir, "cogaps_cancercells.rds")

# Enregistrer l'objet Seurat subset
saveRDS(cogaps_subset, file = output_file)

cat("Objet Seurat filtré et sauvegardé avec succès.\n")
