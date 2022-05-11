# This script performs a lightweight feature QC on all datasets under consideration.

# set directories
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
schraivogel_dir <- .get_config_path("LOCAL_SCHRAIVOGEL_2020_DATA_DIR")
papalexi_dir <- .get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")
liscovitch_dir <- .get_config_path("LOCAL_LISCOVITCH_2021_DATA_DIR")

# load packages
library(ondisc)

# set gene expression threshold
EXPRESSION_THRESH <- 0.005

# Load the datasets
# i. Papalexi gene
papalexi_gene <- read_odm(odm_fp = paste0(papalexi_dir, "processed/gene/expression_matrix.odm"),
                          metadata_fp = paste0(papalexi_dir, "processed/gene/metadata.rds"))

# ii. Schraivogel TAP
schraivogel_tap <- read_odm(odm_fp = paste0(schraivogel_dir, "processed/ground_truth_tapseq/gene/expression_matrix.odm"),
                            metadata_fp = paste0(schraivogel_dir, "processed/ground_truth_tapseq/gene/metadata.rds"))

# iii. Schraivogel Perturb
schraivogel_perturb <- read_odm(odm_fp = paste0(schraivogel_dir, "/processed/ground_truth_perturbseq/gene/expression_matrix.odm"),
                                metadata_fp = paste0(schraivogel_dir, "/processed/ground_truth_perturbseq/gene/metadata.rds"))

# do feature selection Papalexi gene, schraivogel tap, and schraivogel perturb; keep those genes with expression >= 0.005.
filter_features <- function(odm) {
  odm[get_highly_expressed_features(odm, frac_expressed = 0.005),]
}
papalexi_gene_filtered <- filter_features(papalexi_gene)
schraivogel_tap_filtered <- filter_features(schraivogel_tap)
schraivogel_perturb_filtered <- filter_features(schraivogel_perturb)

# save subsetted ODMs
# papalexi gene
save_odm(odm = papalexi_gene_filtered,
         metadata_fp = paste0(sceptre2_dir, "data/papalexi/gene/filtered_metadata.rds"))
# schraivogel tap
save_odm(odm = schraivogel_tap_filtered,
         metadata_fp = paste0(sceptre2_dir, "data/schraivogel/tap/filtered_metadata.rds"))
# schraivogel perturb
save_odm(odm = schraivogel_perturb_filtered,
         metadata_fp = paste0(sceptre2_dir, "data/schraivogel/perturb/filtered_metadata.rds"))
# liscovitch exp 1, atac modality
save_odm(odm = get_modality(exp_1_multimodal_qc, "atac"),
         metadata_fp = paste0(sceptre2_dir, "data/liscovitch/experiment_1/chip_count/filtered_metadata.rds"))
# liscovitch exp 1, gRNA modality
save_odm(odm = get_modality(exp_1_multimodal_qc, "gRNA"),
         metadata_fp = paste0(sceptre2_dir, "data/liscovitch/experiment_1/gRNA/filtered_metadata.rds"))
# liscovitch exp 2, atac modality
save_odm(odm = get_modality(exp_2_multimodal_qc, "atac"),
         metadata_fp = paste0(sceptre2_dir, "data/liscovitch/experiment_2/chip_count/filtered_metadata.rds"))
# liscovitch exp 2, gRNA modality
save_odm(odm = get_modality(exp_2_multimodal_qc, "gRNA"),
         metadata_fp = paste0(sceptre2_dir, "data/liscovitch/experiment_2/gRNA/filtered_metadata.rds"))

