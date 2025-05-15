# Script: CoGAPS Analysis on scRNA-seq Data Using Seurat
#
# Description:
# This script loads a combined Seurat object containing scRNA-seq data,
# extracts the raw count matrix, normalizes it using log1p transformation,
# and performs a CoGAPS (Coordinate Gene Activity in Pattern Sets) analysis
# with parameters optimized for distributed parallel computing.
# The results are saved to a specified output directory.
#
# Note:
# The analysis was run on the Nautilus compute cluster,
# utilizing distributed parallelization to improve computational efficiency.
#
# Inputs:
# - Combined Seurat object RDS file containing raw counts
#
# Outputs:
# - CoGAPS results saved as an RDS file
library(Seurat)
library(SingleCellExperiment)
library(CoGAPS)
library(Matrix)

input_file <- "/scratch/nautilus/users/sanglier-l@univ-nantes.fr/cogaps/data/seurat_combined.rds"
output_dir <- "/scratch/nautilus/users/sanglier-l@univ-nantes.fr/cogaps/output/res_16p"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

pdac_data <- readRDS(input_file)
pdac_epi_counts <- as.matrix(pdac_data@assays$RNA$counts)
norm_pdac_epi_counts <- log1p(pdac_epi_counts)

head(pdac_epi_counts, n = c(5L, 2L))
head(norm_pdac_epi_counts, n = c(5L, 2L))

pdac_params <- CogapsParams(
  nIterations = 50000,  # 50000 itérations
  seed = 42,  # pour la reproductibilité
  nPatterns = 16,  # chaque thread apprendra 20 patterns
  sparseOptimization = TRUE,  # optimiser pour des données éparses
  distributed = "genome-wide"  # parallélisation sur l'ensemble du génome
)

pdac_params
pdac_params <- setDistributedParams(pdac_params, nSets = 7)
pdac_params

pdac_epi_result <- CoGAPS(norm_pdac_epi_counts, pdac_params)

output_file <- file.path(output_dir, "cogaps.rds")
saveRDS(pdac_epi_result, output_file)

