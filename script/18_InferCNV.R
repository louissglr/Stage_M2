library(Seurat)
combine_sample_dir <- "/scratch/nautilus/users/sanglier-l@univ-nantes.fr/infer/data"
rds_path <- file.path(combine_sample_dir, "seurat_cogaps_res_9p.rds")

seurat_combined <- readRDS(rds_path)

library(infercnv)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=GetAssayData(seurat_combined, assay = "RNA", slot = "counts"),
                                    annotations_file=as.matrix(seurat_combined$pattern_cogaps),
                                    delim="\t",
                                    gene_order_file="/scratch/nautilus/users/sanglier-l@univ-nantes.fr/infer/data/mouse_gencode.GRCm38.p6.vM25.basic.annotation.by_gene_name.infercnv_positions.txt",
                                    ref_group_names=c("Pattern_4","Pattern_5"))  

setwd("/scratch/nautilus/users/sanglier-l@univ-nantes.fr/infer")
output_dir <- "n1n2_res0.6"


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.01, 
                             out_dir=output_dir, 
                             num_threads = 10,
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             leiden_resolution=0.6,
                             leiden_method="simple",
                             leiden_function="modularity",
                             no_prelim_plot = TRUE,
                             plot_steps = FALSE,
                             HMM=TRUE) 
