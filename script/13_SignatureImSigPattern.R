library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)

# === Chargement des données ===
seurat_obj <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")
signature_imsig <- read.table("C:/Users/louis/Desktop/Stage/scrnaseq/signatures/ImSig2_mouse.txt", header = FALSE, stringsAsFactors = FALSE)$V1

# === Calcul du score de signature ===
DefaultAssay(seurat_obj) <- "RNA" 
seurat_obj <- AddModuleScore(seurat_obj, features = list(signature_imsig), name = "ImSig_Score")

# === Thème noir personnalisé ===
theme_black_background <- theme_minimal() +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    axis.line = element_line(color = "white"),
    axis.ticks = element_line(color = "white"),
    axis.text = element_text(color = "white", size = 14),
    axis.title = element_text(color = "white"),
    panel.grid = element_blank(),
    legend.background = element_rect(fill = "black", color = NA),
    legend.box.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    plot.title = element_text(color = "white", hjust = 0.5)
  )

# === Fonction pour ajouter le nombre de cellules ===
add_cell_count <- function(p, data, group_by) {
  counts <- data %>%
    group_by(.data[[group_by]]) %>%
    summarise(
      count = n(),
      max_val = max(ImSig_Score1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_pos = max_val + 0.05)  # marge au-dessus du max
  
  p + geom_text(
    data = counts,
    aes(x = .data[[group_by]], y = y_pos, label = paste("n =", count)),
    color = "white",
    size = 2.7,
    fontface = "bold",
    inherit.aes = FALSE
  )
}

# === Violin plot par pattern CoGAPS ===
p1 <- VlnPlot(
  object = seurat_obj,
  features = "ImSig_Score1",
  group.by = "pattern_cogaps",
  pt.size = 0
) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
  labs(title = "ImSig Score by CoGAPS patterns") +
  theme_black_background +
  theme(legend.position = "none")

# === Calcul de la valeur max pour l’axe Y ===
max_y1 <- max(seurat_obj@meta.data$ImSig_Score1, na.rm = TRUE)

# === Ajustement des axes et étiquettes ===
p1 <- p1 +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(NA, max_y1 * 1.15)) +
  scale_x_discrete(labels = function(x) gsub("Pattern_", "", x))  # Supprime "Pattern_" des labels

# === Ajout du nombre de cellules ===
p1 <- add_cell_count(p1, seurat_obj@meta.data, "pattern_cogaps")

# === Sauvegarde du graphique ===
ggsave(
  filename = "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/13_SignatureImSigPattern/ImSig_ViolinPlot_CoGAPSPatterns.png",
  plot = p1,
  width = 6,
  height = 6,
  dpi = 300,
  bg = "black"
)
