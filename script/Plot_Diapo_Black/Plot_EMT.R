############### Chargement des packages #########################################
library(Seurat)
library(CoGAPS)
library(dplyr)
library(ggplot2)
library(forcats)
library(purrr)
library(tidyr)
library(gtools)  # Pour mixedsort()
library(viridis)

############### Chargement des fichiers Hallmark ###############################
hallmark_path <- "C:/Users/louis/Desktop/Stage/scrnaseq/hallmark"
files <- list.files(hallmark_path, pattern = "\\.txt$", full.names = TRUE)

hallmark_ls <- lapply(files, function(file) {
  genes_line <- readLines(file)
  genes <- unlist(strsplit(genes_line, ","))
  genes <- gsub("\"", "", genes)  # Nettoyage
  genes
})
names(hallmark_ls) <- tools::file_path_sans_ext(basename(files))

############### Chargement objet Seurat pour récupérer les couleurs #############
seurat_for_colors <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/09_AddCoGAPSToSeurat/seurat_cogaps_res_9p.rds")

# Nettoyer noms patterns et définir ordre décroissant (9 → 1)
seurat_for_colors$pattern_cogaps <- gsub("Pattern_", "", seurat_for_colors$pattern_cogaps)
pattern_levels_dec <- paste0("Pattern_", sort(unique(seurat_for_colors$pattern_cogaps), decreasing = TRUE))

# Palette de couleurs UMAP pour les patterns, nommée avec les patterns dans l’ordre décroissant
umap_colors <- scales::hue_pal()(length(pattern_levels_dec))
names(umap_colors) <- pattern_levels_dec

############### Parcours des objets CoGAPS ######################################
res_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/08_CoGAPS/"
subdirs <- mixedsort(list.dirs(res_dir, recursive = FALSE))

seurat_combined <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/04_MergedData/seurat_combined.rds")

for (subdir in subdirs) {
  sample_name <- basename(subdir)
  
  # On ne travaille que sur "res_9p"
  if (sample_name != "res_9p") next
  
  cogaps_file <- file.path(subdir, "cogaps.rds")
  
  if (file.exists(cogaps_file)) {
    message("▶️ Traitement de : ", cogaps_file)
    
    cogaps <- readRDS(cogaps_file)
    seurat <- seurat_combined
    
    # Ajout des patterns comme assay
    patterns_in_order <- t(cogaps@sampleFactors[colnames(seurat), , drop = FALSE])
    seurat[["CoGAPS"]] <- CreateAssayObject(counts = patterns_in_order)
    DefaultAssay(seurat) <- "CoGAPS"
    
    # Analyse ORA
    hallmarks_ora <- getPatternGeneSet(cogaps, gene.sets = hallmark_ls, method = "overrepresentation")
    
    hallmarks_ora_annotated <- lapply(seq_along(hallmarks_ora), function(i) {
      df <- hallmarks_ora[[i]]
      df$Pattern <- paste0("Pattern_", i)
      return(df)
    })
    
    hallmarks_df <- bind_rows(hallmarks_ora_annotated) %>%
      mutate(minus_log10_padj = -10 * log10(padj))
    
    # Barplot ciblé EMT pour patterns 2,3,7,8,9
    patterns_cibles <- paste0("Pattern_", c(2, 3, 7, 8, 9))
    
    # Ordre décroissant dans le plot: 9,8,7,3,2 (par exemple)
    patterns_cibles_ordered <- rev(patterns_cibles)  # Donc ici 9,8,7,3,2
    
    signature_target <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
    
    emt_df <- hallmarks_df %>%
      filter(Pattern %in% patterns_cibles) %>%
      filter(gene.set == signature_target) %>%
      mutate(Pattern = factor(Pattern, levels = patterns_cibles_ordered))
    
    # Couleurs selon palette Seurat pour patterns_cibles_ordered
    pattern_colors_to_use <- umap_colors[patterns_cibles_ordered]
    
    emt_plot <- ggplot(emt_df, aes(x = minus_log10_padj, y = Pattern, fill = Pattern)) +
      geom_col(width = 0.6) +
      geom_vline(xintercept = -10 * log10(0.05), linetype = "dotted", color = "red") +
      geom_text(aes(label = signif(padj, 2), color = padj < 0.05), hjust = -0.2, size = 6) +
      scale_color_manual(values = c("white", "red")) +
      scale_fill_manual(values = pattern_colors_to_use) +
      theme_minimal(base_family = "sans") +
      theme(
        panel.background = element_rect(fill = "black", color = NA),
        plot.background = element_rect(fill = "black", color = NA),
        legend.background = element_rect(fill = "black", color = NA),
        legend.box.background = element_rect(fill = "black", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "white", size = 14),
        axis.title = element_text(color = "white", size = 16, face = "bold"),
        legend.text = element_text(color = "white", size = 14),
        legend.title = element_text(color = "white", size = 16, face = "bold"),
        axis.line = element_line(color = "white"),
        axis.ticks = element_line(color = "white"),
        plot.title = element_text(hjust = 0.5, color = "white", size = 18, face = "bold"),
        legend.position = "none"
      ) +
      labs(
        title = paste0("EMT Enrichment - ", sample_name),
        x = expression(-10*log[10](padj)),
        y = "Pattern"
      ) +
      coord_cartesian(xlim = c(0, max(emt_df$minus_log10_padj, na.rm = TRUE) + 15))
    
    print(emt_plot)
    
  } else {
    warning("Fichier cogaps.rds manquant dans : ", subdir)
  }
}

ggsave("C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Diapo_Black/EMT_Plot.png",
       emt_plot, width = 18, height = 6, dpi = 300, bg = "black")
