# This script uses sceptre3 to run an approximate undercover analysis on the 
LOCAL_SCEPTRE2_DATA_DIR <-.get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
result_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "results/writeup_results/calibration_check")
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
ifn_gamma_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "data/frangieh/ifn_gamma/")

# gene info
gene_odm_fp <- paste0(ifn_gamma_dir, "gene/matrix.odm")

# grna info
grna_odm_fp <- paste0(ifn_gamma_dir, "grna_assignment/matrix.odm")

# mm odm metadata fp
mm_metadata_fp <- paste0(ifn_gamma_dir, "multimodal_metadata.rds")

# construct mm odm
mm_odm <- ondisc::read_multimodal_odm(odm_fps = c(gene_odm_fp, grna_odm_fp),
                                      multimodal_metadata_fp = mm_metadata_fp)

# get the in-memory gene matrix
gene_odm <- mm_odm |> ondisc::get_modality("gene")
response_matrix <- gene_odm[[seq(1, nrow(gene_odm)),]]
rownames(response_matrix) <- ondisc::get_feature_ids(gene_odm)

# get the in-memory grna matrix
grna_odm <- mm_odm |> ondisc::get_modality("grna_assignment")
grna_matrix <- grna_odm[[seq(1, nrow(grna_odm)),]]

# covariate matrix
covariate_data_frame <- mm_odm |> ondisc::get_cell_covariates()

# grna group data frame
grna_group_data_frame <- data.frame(grna_id = rownames(grna_odm@feature_covariates),
                                    grna_group = grna_odm@feature_covariates$target)
n_ntc <- grna_group_data_frame |>
  dplyr::filter(grna_group == "non-targeting") |>
  nrow()

# set formulas, grna group target name
gene_formula <- mm_odm@global_misc$formula

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
                                                  n_calibration_pairs = 200) # n_ntc * nrow(response_matrix))

saveRDS(object = undercover_result, file = paste0(result_dir, "/frangieh_calibration_check.rds"))
