# load libraries and resolve conflicts
library(ondisc) 
library(sceptre)
library(readr)
library(dplyr)
library(conflicted)
conflicted::conflicts_prefer(dplyr::filter)

# set up directories
LOCAL_SCEPTRE2_DATA_DIR <-.get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
frangieh_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, 
                               "data/frangieh/control/")

# gene info
gene_odm_fp <- paste0(frangieh_dir, "gene/matrix.odm")
gene_metadata_fp <- paste0(frangieh_dir, "gene/metadata_qc.rds")
gene_odm <- read_odm(odm_fp = gene_odm_fp, metadata_fp = gene_metadata_fp)
gene_covariate_matrix <- gene_odm |> get_cell_covariates() 
gene_expression_matrix <- gene_odm[[seq(1, nrow(gene_odm)),]]
rownames(gene_expression_matrix) <- get_feature_ids(gene_odm)

# grna info
grna_odm_fp <- paste0(frangieh_dir, "grna_assignment/matrix.odm")
grna_metadata_fp <- paste0(frangieh_dir, "grna_assignment/metadata_qc.rds")
grna_odm <- read_odm(odm_fp = grna_odm_fp, metadata_fp = grna_metadata_fp)
grna_matrix <- grna_odm[[seq(1, nrow(grna_odm)),]]
grna_groups <- data.frame(grna_id = rownames(grna_odm@feature_covariates),
                          grna_group = grna_odm@feature_covariates$target)

# set arguments for SCEPTRE
response_matrix <- gene_expression_matrix
grna_matrix <- grna_matrix
rownames(grna_matrix) <- ondisc::get_feature_ids(grna_odm)
covariate_data_frame <- gene_covariate_matrix
grna_group_data_frame <- grna_groups
formula_object <- ~log(n_umis) + log(n_nonzero)
calibration_check <- FALSE
unique_grna <- unique(grna_groups$grna_group)
response_grna_group_pairs <- expand.grid(response_id = get_feature_ids(gene_odm),
                                         grna_group = unique_grna[-which(unique_grna == 'non-targeting')])

# run SCEPTRE
result_sceptre <- run_sceptre_lowmoi(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  covariate_data_frame = covariate_data_frame,
  grna_group_data_frame = grna_group_data_frame,
  formula_object = formula_object,
  response_grna_group_pairs = response_grna_group_pairs,
  calibration_check = calibration_check,
  return_debugging_metrics = TRUE
)

# save the results
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
output_dir <- paste0(sceptre2_dir, "results/frangieh_analysis/")
if (!dir.exists(output_dir)) dir.create(output_dir)
saveRDS(result_sceptre, paste0(output_dir, "sceptre_frangieh_control_results.rds"))