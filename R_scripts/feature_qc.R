# set directories
# SCEPTRE 2 DIR
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
# SCHRAIVOGEL DIR
schraivogel_dir <- .get_config_path("LOCAL_SCHRAIVOGEL_2020_DATA_DIR")
# PAPALEXI DIR
papalexi_dir <- .get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")

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

# do feature selection on all datasets; keep those genes with expression >= 0.005.
filter_features <- function(odm) {
  odm[get_highly_expressed_features(odm, frac_expressed = 0.005),]
}
papalexi_gene_filtered <- filter_features(papalexi_gene)
schraivogel_tap_filtered <- filter_features(schraivogel_tap)
schraivogel_perturb_filtered <- filter_features(schraivogel_perturb)

# save subsetted ODMs
save_odm(odm = papalexi_gene_filtered,
         metadata_fp = paste0(sceptre2_dir, "data/papalexi/gene/filtered_metadata.rds"))
save_odm(odm = schraivogel_tap_filtered,
         metadata_fp = paste0(sceptre2_dir, "data/schraivogel/tap/filtered_metadata.rds"))
save_odm(odm = schraivogel_perturb_filtered,
         metadata_fp = paste0(sceptre2_dir, "data/schraivogel/perturb/filtered_metadata.rds"))
