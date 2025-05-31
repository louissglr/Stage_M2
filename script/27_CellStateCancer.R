############### Packages ########################################################
library(Seurat)
library(CoGAPS)
library(dplyr)
library(ggplot2)
library(forcats)
library(purrr)
library(tidyr)
library(gtools)
library(tools)
library(tidytext)   # pour reorder_within() et scale_y_reordered()
library(openxlsx)   # pour l'export Excel
library(ggtext)     # pour colorer les titres des facettes

############### Chargement des Hallmarks #######################################
hallmark_path <- "C:/Users/louis/Desktop/Stage/scrnaseq/hallmark"
files <- list.files(hallmark_path, pattern = "\\.txt$", full.names = TRUE)

hallmark_ls <- lapply(files, function(file) {
  genes_line <- readLines(file)
  genes <- unlist(strsplit(genes_line, ","))
  genes <- gsub("\"", "", genes)
  genes
})
names(hallmark_ls) <- file_path_sans_ext(basename(files))

############### Chargement Seurat & CoGAPS #####################################
seurat_combined <- readRDS("C:/Users/louis/Desktop/Stage/scrnaseq_final/output/19_Keep_CoGAPS_Cancer/cogaps_cancercells.rds")

res_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/20_CoGAPS_cancer/"
output_dir <- "C:/Users/louis/Desktop/Stage/scrnaseq_final/output/24_HallmarkTop3_Figure/"

subdirs <- mixedsort(list.dirs(res_dir, recursive = FALSE))

for (subdir in subdirs) {
  cogaps_file <- file.path(subdir, "cogaps.rds")
  
  if (file.exists(cogaps_file)) {
    message("▶️ Traitement de : ", cogaps_file)
    
    sample_name <- basename(subdir)
    cogaps <- readRDS(cogaps_file)
    seurat <- seurat_combined
    
    # Ajouter les patterns à l'objet Seurat
    patterns <- t(cogaps@sampleFactors[colnames(seurat), , drop = FALSE])
    seurat[["CoGAPS"]] <- CreateAssayObject(counts = patterns)
    DefaultAssay(seurat) <- "CoGAPS"
    
    # ORA Hallmark
    hallmarks_ora <- getPatternGeneSet(cogaps, gene.sets = hallmark_ls, method = "overrepresentation")
    
    hallmarks_df <- map2_dfr(hallmarks_ora, seq_along(hallmarks_ora), ~ {
      .x$Pattern <- paste0("Pattern_", .y)
      .x
    }) %>%
      mutate(minus_log10_padj = -10 * log10(padj))
    
    # Sélection des top 3 enrichissements par pattern
    top3_df <- hallmarks_df %>%
      group_by(Pattern) %>%
      arrange(padj) %>%
      slice_head(n = 3) %>%
      ungroup()
    
    # Définir les couleurs pour chaque Pattern
    pattern_colors <- c(
      "Pattern_1" = "#f465e3",
      "Pattern_2" = "#609dff",
      "Pattern_3" = "#01bec4",
      "Pattern_4" = "#00ba38",
      "Pattern_5" = "#b79f00",
      "Pattern_6" = "#f8776c"
    )
    
    # Ajouter des balises HTML pour colorer les titres
    top3_df <- top3_df %>%
      mutate(Pattern_colored = paste0(
        "<span style='color:", pattern_colors[Pattern], "'>", Pattern, "</span>"
      ))
    
    top3_df$Pattern_colored <- factor(top3_df$Pattern_colored, levels = unique(top3_df$Pattern_colored))
    
    # Figure : facettes ggplot avec ordre local par pattern
    p <- ggplot(top3_df, aes(x = minus_log10_padj,
                             y = reorder_within(gene.set, -padj, Pattern_colored),
                             fill = minus_log10_padj)) +
      geom_col(width = 0.6) +
      facet_wrap(~ Pattern_colored, scales = "free_y", ncol = 2) +
      scale_fill_viridis_c(option = "C") +
      scale_y_reordered() +
      theme_minimal(base_size = 13) +
      xlab(expression(-10*log[10](padj))) +
      ylab("Hallmark") +
      geom_vline(xintercept = -10 * log10(0.05), linetype = "dotted", color = "red") +
      geom_text(aes(label = signif(padj, 2), color = padj < 0.05),
                hjust = -0.2, size = 3) +
      scale_color_manual(values = c("black", "red")) +
      theme(
        legend.position = "none",
        strip.text = ggtext::element_markdown(
          size = 14,
          face = "bold",
          margin = margin(5, 5, 5, 5)
        )
        
        ,
        plot.margin = margin(10, 10, 10, 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")
      ) +
      xlim(0, max(top3_df$minus_log10_padj, na.rm = TRUE) + 10)
    
    # Sauvegarde PNG
    output_file <- file.path(output_dir, paste0("hallmark_top3_", sample_name, ".png"))
    if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
    
    ggsave(output_file, plot = p, width = 14, height = max(6, 2.5 * length(unique(top3_df$Pattern)) / 2),
           dpi = 300, bg = "white")
    
    message("✅ Image enregistrée : ", output_file)
    
    # Export Excel si le dossier est "res_6p"
    if (basename(subdir) == "res_6p") {
      message("📁 Export Excel pour le répertoire : res_6p")
      excel_df <- top3_df[, c("gene.set", "pval", "padj", "Pattern")]
      excel_output_file <- file.path(output_dir, "res_6patterns_res6p.xlsx")
      write.xlsx(excel_df, file = excel_output_file, sheetName = "Top3_Hallmarks", overwrite = TRUE)
    }
    
  } else {
    warning("❌ Fichier manquant : ", cogaps_file)
  }
}
