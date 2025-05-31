# Charger l'objet Seurat
seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq/output/03_CombineSample/seurat_combined.rds")

# Charger les données d'expression normalisée
expr <- GetAssayData(seurat_obj, slot = "data")

# Nom exact des gènes (adapter si différent)
gene_td <- "TdTomato"
gene_sv <- "SV40"

# Identifier les cellules positives (> 0) pour chaque gène
TdTomato_pos <- colnames(seurat_obj)[expr[gene_td, ] > 0]
SV40_pos     <- colnames(seurat_obj)[expr[gene_sv, ] > 0]

# Extraire les métadonnées
meta <- seurat_obj@meta.data

# Ajouter une colonne identifiant TdTomato et SV40 positif
meta$TdTomato_pos <- rownames(meta) %in% TdTomato_pos
meta$SV40_pos     <- rownames(meta) %in% SV40_pos

# Compter les cellules positives par orig.ident
library(dplyr)

count_table <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    TdTomato = sum(TdTomato_pos),
    SV40 = sum(SV40_pos)
  ) %>%
  t()  # transpose

# Donner des noms clairs
colnames(count_table) <- count_table[1, ]
count_table <- count_table[-1, ]

# Afficher le tableau final
print(count_table)

# Charger dplyr
library(dplyr)

# Compter le total de cellules par échantillon
total_cells <- meta %>%
  group_by(orig.ident) %>%
  summarise(total = n()) %>%
  as.data.frame()

# Compter les cellules positives à TdTomato et SV40
positive_counts <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    TdTomato = sum(TdTomato_pos),
    SV40 = sum(SV40_pos)
  ) %>%
  as.data.frame()

# Ajouter le total dans le même tableau
merged <- left_join(positive_counts, total_cells, by = "orig.ident")

# Ajouter les pourcentages au format "count (xx%)"
format_pct <- function(count, total) {
  pct <- round(100 * count / total, 1)
  paste0(count, " (", pct, "%)")
}

# Appliquer la fonction à chaque colonne de marqueurs
formatted <- data.frame(
  Marker = c("TdTomato", "SV40"),
  t(sapply(c("TdTomato", "SV40"), function(marker) {
    mapply(format_pct, merged[[marker]], merged$total)
  }))
)

# Mettre les noms des colonnes correctement
colnames(formatted) <- c("Marker", merged$orig.ident)

# Afficher le tableau final
print(formatted, row.names = FALSE)

