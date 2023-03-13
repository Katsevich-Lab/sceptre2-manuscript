
library(ondisc) # devtools::install_github('timothy-barry/ondisc')
library(sceptre3)
library(BH)
# This script uses sceptre3 to run an approximate undercover analysis on the 
LOCAL_SCEPTRE2_DATA_DIR <-.get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
papalexi_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "data/papalexi/eccite_screen/")

# gene info
gene_odm_fp <- paste0(papalexi_dir, "gene/matrix.odm")

# grna info
grna_odm_fp <- paste0(papalexi_dir, "grna_expression/matrix.odm")

# mm odm metadata fp
mm_metadata_fp <- paste0(papalexi_dir, "multimodal_metadata.rds")

# construct mm odm
mm_odm <- ondisc::read_multimodal_odm(odm_fps = c(gene_odm_fp, grna_odm_fp),
                                      multimodal_metadata_fp = mm_metadata_fp)

# get the in-memory gene matrix
gene_odm <- mm_odm |> ondisc::get_modality("gene")
response_matrix <- gene_odm[[seq(1, nrow(gene_odm)),]]
rownames(response_matrix) <- ondisc::get_feature_ids(gene_odm)

# get the in-memory grna matrix
grna_odm <- mm_odm |> ondisc::get_modality("grna_expression")
grna_matrix <- grna_odm[[seq(1, nrow(grna_odm)),]]

# covariate matrix
covariate_data_frame <- mm_odm |> ondisc::get_cell_covariates()

# grna group data frame
grna_group_data_frame <- data.frame(grna_id = rownames(grna_odm@feature_covariates),
                                    grna_group = grna_odm@feature_covariates$target)

# set formulas, grna group target name
gene_formula <- ~ log(gene_n_umis) + log(gene_n_nonzero) + bio_rep + phase + p_mito

undercover_result <- sceptre3::run_sceptre_lowmoi(response_matrix = response_matrix,
                                                  grna_matrix = grna_matrix,
                                                  covariate_data_frame = covariate_data_frame,
                                                  grna_group_data_frame = grna_group_data_frame,
                                                  formula_object = gene_formula,
                                                  calibration_check = TRUE,
                                                  response_grna_group_pairs = NULL,
                                                  test_stat = "full",
                                                  return_resampling_dist = FALSE,
                                                  fit_skew_normal = TRUE,
                                                  B1 = 500,
                                                  B2 = 5000,
                                                  B3 = 25000,
                                                  undercover_group_size = 1,
                                                  n_calibration_pairs = 9 * nrow(response_matrix))
p <- sceptre3::plot_calibration_results(undercover_result)
library(ggplot2)
