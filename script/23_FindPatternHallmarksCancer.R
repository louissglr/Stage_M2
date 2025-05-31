# Script: Hallmark Pathway Enrichment per CoGAPS Pattern (ORA + UMAP Visualization)
#
# Description:
# This script performs Hallmark gene set overrepresentation analysis (ORA) on CoGAPS-derived patterns
# for single-cell RNA-seq data and visualizes the enrichment results per pattern, along with CoGAPS
# pattern expression across cells via UMAP. It outputs a PDF per sample.
#
# Inputs:
# - Hallmark gene sets as `.txt` files (comma-separated gene names).
# - CoGAPS results (`cogaps.rds`) from subdirectories inside a main results folder.
# - A Seurat object (`seurat_combined.rds`) with UMAP embeddings.
#
# Outputs:
# - A PDF file per sample containing:
#     • FeaturePlots for each CoGAPS pattern on UMAP
#     • Barplots of Hallmark pathway enrichment per pattern
#
# Key Steps:
# - Load and clean Hallmark gene sets from plain-text files.
# - For each CoGAPS result:
#     • Extract pattern weights and add them to the Seurat object as a new assay.
#     • Visualize each pattern on UMAP using `FeaturePlot`.
#     • Perform ORA using `getPatternGeneSet` from CoGAPS.
#     • Create annotated barplots of significant enrichments.
#     • Save results to a PDF
############### Chargement des packages #########################################
library(Seurat)
library(CoGAPS)
library(dplyr)
library(ggplot2)
library(forcats)
library(purrr)
library(tidyr)
library(gtools)  # Pour mixedsort()

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

############### Parcours des objets CoGAPS ######################################
res_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/20_CoGAPS_cancer/"
output_base <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/23_FindPatternHallmarksCancer/"

subdirs <- mixedsort(list.dirs(res_dir, recursive = FALSE))

seurat_combined <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/19_Keep_CoGAPS_Cancer/cogaps_cancercells.rds")

for (subdir in subdirs) {
  cogaps_file <- file.path(subdir, "cogaps.rds")
  
  if (file.exists(cogaps_file)) {
    message("▶️ Traitement de : ", cogaps_file)
    
    sample_name <- basename(subdir)
    
    cogaps <- readRDS(cogaps_file)
    
    seurat <- seurat_combined
    patterns_in_order <- t(cogaps@sampleFactors[colnames(seurat), , drop = FALSE])
    seurat[["CoGAPS"]] <- CreateAssayObject(counts = patterns_in_order)
    DefaultAssay(seurat) <- "CoGAPS"
    pattern_names <- rownames(seurat[["CoGAPS"]])
    
    featureplot_obj <- FeaturePlot(seurat, pattern_names, reduction = "umap")
    
    
    # Analyse ORA
    hallmarks_ora <- getPatternGeneSet(cogaps, gene.sets = hallmark_ls, method = "overrepresentation")
    
    hallmarks_ora_annotated <- lapply(seq_along(hallmarks_ora), function(i) {
      df <- hallmarks_ora[[i]]
      df$Pattern <- paste0("Pattern_", i)
      return(df)
    })
    
    hallmarks_df <- bind_rows(hallmarks_ora_annotated) %>%
      mutate(minus_log10_padj = -10 * log10(padj))
    
    pattern_plot_list <- hallmarks_df %>%
      split(.$Pattern) %>%
      map(~ {
        df <- .x %>%
          mutate(gene.set = fct_reorder(gene.set, minus_log10_padj))
        
        ggplot(df, aes(x = minus_log10_padj, y = gene.set, fill = minus_log10_padj)) +
          geom_col() +
          scale_fill_viridis_c(option = "C") +
          theme_minimal() +
          ggtitle(paste0("Hallmark Enrichment - ", unique(df$Pattern))) +
          xlab(expression(-10*log[10](padj))) +
          ylab("Hallmark") +
          geom_vline(xintercept = -10 * log10(0.05), linetype = "dotted", color = "red") +
          geom_text(aes(label = signif(padj, 2), color = padj < 0.05), hjust = -0.2, size = 3) +
          scale_color_manual(values = c("black", "red")) +
          theme(legend.position = "none") +
          xlim(0, max(df$minus_log10_padj, na.rm = TRUE) + 10)
      })
    
    output_pdf <- file.path(output_base, paste0("hallmarks_", sample_name, ".pdf"))
    if (!dir.exists(dirname(output_pdf))) {
      dir.create(dirname(output_pdf), recursive = TRUE)
    }
    
    pdf(output_pdf, width = 10, height = 7)
    print(featureplot_obj)  
    walk(pattern_plot_list, print)  
    dev.off()
    
    message("PDF généré : ", output_pdf)
    
  } else {
    warning("Fichier cogaps.rds manquant dans : ", subdir)
  }
}
