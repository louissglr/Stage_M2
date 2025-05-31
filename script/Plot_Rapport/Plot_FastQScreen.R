library(tidyverse)
library(stringr)
library(purrr)

read_fastq_screen <- function(file) {
  cat("Lecture de :", file, "\n")
  
  # Lire toutes les lignes du fichier
  lines <- readLines(file)
  
  # Extraire la valeur %Hit_no_genomes
  hit_no_genomes_line <- grep("^%Hit_no_genomes:", lines, value = TRUE)
  no_hits_value <- if(length(hit_no_genomes_line) == 1) {
    as.numeric(str_extract(hit_no_genomes_line, "\\d+\\.\\d+"))
  } else {
    NA_real_
  }
  
  # Trouver l'index de la ligne d'en-tête du tableau (celle qui commence par "Genome")
  header_index <- grep("^Genome", lines)
  
  # Trouver l'index de la ligne %Hit_no_genomes (fin du tableau)
  hit_index <- grep("^%Hit_no_genomes:", lines)
  
  # Extraire uniquement les lignes du tableau, sans la ligne %Hit_no_genomes
  table_lines <- lines[header_index:(hit_index - 1)]
  
  # Lire le tableau depuis ces lignes
  df <- read.delim(text = table_lines, check.names = FALSE)
  
  required_cols <- c(
    "Genome",
    "%One_hit_one_genome",
    "%Multiple_hits_one_genome",
    "%One_hit_multiple_genomes",
    "%Multiple_hits_multiple_genomes"
  )
  
  missing_cols <- setdiff(required_cols, colnames(df))
  if(length(missing_cols) > 0) {
    warning("Colonnes manquantes dans : ", file, "\nManquantes : ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }
  
  # Extraction du nom brut de l’échantillon
  Sample_raw <- basename(file) %>% str_remove("_R2_001_screen.txt")
  
  # Renommage des échantillons selon ta demande
  Sample <- case_when(
    Sample_raw == "N1-0h-2-3"    ~ "N2-0h",
    Sample_raw == "N1-2week-1"   ~ "N1-2week",
    Sample_raw == "N2-1week-1-2" ~ "N2-1week",
    Sample_raw == "N1-0h-1"      ~ "N1-0h",
    Sample_raw == "N1-4week-1-2" ~ "N1-4week",
    Sample_raw == "N2-72h-1-2"   ~ "N2-72h",
    Sample_raw == "N1-1week-1-2" ~ "N1-1week",
    Sample_raw == "N1-72h-1"     ~ "N1-72h",
    TRUE                         ~ Sample_raw
  )
  
  df_long <- df %>%
    select(all_of(required_cols)) %>%
    pivot_longer(
      cols = -Genome,
      names_to = "Hit_Type",
      values_to = "Percentage"
    ) %>%
    mutate(
      Sample = Sample,
      Hit_Type = case_when(
        Hit_Type == "%One_hit_one_genome" ~ "One hit / one genome",
        Hit_Type == "%Multiple_hits_one_genome" ~ "Multiple hits / one genome",
        Hit_Type == "%One_hit_multiple_genomes" ~ "One hit / multiple genomes",
        Hit_Type == "%Multiple_hits_multiple_genomes" ~ "Multiple hits / multiple genomes",
        TRUE ~ Hit_Type
      )
    )
  
  # Ajouter une ligne "No hits" par échantillon avec la valeur extraite
  no_hits_df <- tibble(
    Genome = "No hits",
    Hit_Type = "No hits",
    Percentage = no_hits_value,
    Sample = Sample
  )
  
  bind_rows(df_long, no_hits_df)
}

# Chemin vers le dossier contenant les fichiers
folder_path <- "C:/Users/louis/Desktop/Stage/scrnaseq/fastq_screen_results"

# Liste des fichiers
files <- list.files(folder_path, pattern = "_screen.txt$", full.names = TRUE)

genome_levels <- c(
  "No hits", "Adapters", "Vectors", "Lambda", "PhiX", "MT", "rRNA",
  "Ecoli", "Arabidopsis", "Yeast", "Worm", "Drosophila", "Rat",
  "Human", "Mouse"
)

sample_levels <- c(
  "N1-0h", "N1-72h", "N1-1week", "N1-2week", "N1-4week",
  "N2-0h", "N2-72h", "N2-1week"
)

all_data <- all_data %>%
  mutate(Sample = factor(Sample, levels = sample_levels))



# Plot barplot empilé horizontal, proportions sur 0-100, panel par sample
p <- ggplot(all_data, aes(x = Genome, y = Percentage, fill = Hit_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  facet_wrap(~ Sample, scales = "free_y") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_fill_manual(
    values = c(
      "No hits" = "grey",
      "Multiple hits / multiple genomes" = "red",
      "One hit / multiple genomes" = "orange",
      "Multiple hits / one genome" = "blue",
      "One hit / one genome" = "green"
    )
  ) +
  labs(
    title = "FastQ Screen hits proportions",
    y = "Percentage (%)",
    x = "Genome",
    fill = "Hit type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom"
  )
p <- p + theme(
  plot.title = element_text(size = 20, face = "bold"),        # Titre principal plus gros et en gras
  axis.title.x = element_text(size = 16),                     # Label axe x plus gros
  axis.title.y = element_text(size = 16),                     # Label axe y plus gros
  strip.text = element_text(size = 14, face = "bold"),        # Titre des facets (Sample) plus gros et gras
  legend.title = element_text(size = 14),                     # Titre légende plus gros
  legend.text = element_text(size = 12)                       # Texte légende un peu plus grand
)

ggsave(
  filename = "fastq_screen_hits.png",
  plot = p,
  path = "C:/Users/louis/Desktop/Stage/scrnaseq_final/script/Plot_Rapport",
  width = 12, height = 8, units = "in", dpi = 300,
  bg = "white"
)

