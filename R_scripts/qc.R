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
# Liscovitch QC function
run_lisc_qc <- function(multimodal_odm, n_gRNA_reads = 100, n_atac_frags = 500, p_reads_one_gRNA_seq) {
  gRNA_mod <- get_modality(multimodal_odm, "gRNA")
  assignable <- gRNA_mod[[1:nrow(gRNA_mod),]] |> apply(MARGIN = 2, FUN = function(col) {
    max(col)/sum(col) > p_reads_one_gRNA_seq
  })
  multimodal_odm <- multimodal_odm[,assignable]
  high_quality <-  multimodal_odm |>
    get_cell_covariates() |> 
    dplyr::filter(atac_n_fragments >= n_atac_frags) |>
    dplyr::filter(gRNA_n_umis >= 100) |> row.names()
  multimodal_odm[,high_quality]
}

# iv. Liscovitch exp 1 (QC metrics: (a) >= 100 gRNA reads, (b) >= 99% of reads assigned to one gRNA sequence, (c) >= 500 ATAC-seq fragments)
liscovitch_atac_1 <- read_odm(odm_fp = paste0(liscovitch_dir, "processed/experiment_1/chip_counts/chip_counts.odm"),
                            metadata_fp = paste0(liscovitch_dir, "processed/experiment_1/chip_counts/metadata.rds"))
liscovitch_gRNA_1 <- read_odm(odm_fp = paste0(liscovitch_dir, "processed/experiment_1/gRNA_counts/gRNA_counts.odm"),
                            metadata_fp = paste0(liscovitch_dir, "processed/experiment_1/gRNA_counts/metadata.rds"))
exp_1_multimodal <- multimodal_ondisc_matrix(covariate_ondisc_matrix_list = list(atac = liscovitch_atac_1, gRNA = liscovitch_gRNA_1))
exp_1_multimodal_qc <- run_lisc_qc(exp_1_multimodal, p_reads_one_gRNA_seq = 0.99)

# v. Liscovitch exp 2 (QC metrics: >= 100 gRNA reads, (b) >= 90% of reads assigned to one gRNA sequence, (c) >= 500 ATAC-seq fragments)
liscovitch_atac_2 <- read_odm(odm_fp = paste0(liscovitch_dir, "processed/experiment_2/chip_counts/chip_counts.odm"),
                              metadata_fp = paste0(liscovitch_dir, "processed/experiment_2/chip_counts/metadata.rds"))
liscovitch_gRNA_2 <- read_odm(odm_fp = paste0(liscovitch_dir, "processed/experiment_2/gRNA_counts/gRNA_counts.odm"),
                              metadata_fp = paste0(liscovitch_dir, "processed/experiment_2/gRNA_counts/metadata.rds"))
exp_2_multimodal <- multimodal_ondisc_matrix(covariate_ondisc_matrix_list = list(atac = liscovitch_atac_2, gRNA = liscovitch_gRNA_2))
exp_2_multimodal_qc <- run_lisc_qc(exp_2_multimodal, p_reads_one_gRNA_seq = 0.90)


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
save_odm(odm = get_modality(exp_1_multimodal_qc, "atac"),
         metadata_fp = paste0(sceptre2_dir, "data/liscovitch/experiment_1/chip_count/filtered_metadata.rds"))
# liscovitch exp 2, gRNA modality
save_odm(odm = get_modality(exp_1_multimodal_qc, "gRNA"),
         metadata_fp = paste0(sceptre2_dir, "data/liscovitch/experiment_1/gRNA/filtered_metadata.rds"))

