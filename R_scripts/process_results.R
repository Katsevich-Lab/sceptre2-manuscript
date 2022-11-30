library(lowmoi)
library(tidyverse)
sceptre2_results_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
sample_size_df <- readRDS(paste0(sceptre2_results_dir, "dataset_sample_sizes/n_nonzero_cells_per_grna.rds"))

# undercover res grp = 1
undercover_res <- readRDS(paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1.rds"))
undercover_res_extra <- readRDS(paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1_debug.rds")) |>
  dplyr::mutate(method = forcats::fct_recode(method, "sceptre_no_covariates" = "sceptre"))
undercover_res <- rbind(undercover_res, undercover_res_extra)
undercover_res_processed <- process_undercover_result(undercover_res, sample_size_df) |>
  mutate(p_value = ifelse(p_value <= 0, 1e-8, p_value))
saveRDS(undercover_res_processed, paste0(sceptre2_results_dir,
                                         "undercover_grna_analysis/undercover_result_grp_1_processed.rds"))
rm(list = c("undercover_res", "undercover_res_extra", "undercover_res_processed"))

# resampling results
resampling_res <- readRDS(paste0(sceptre2_results_dir,
                                 "resampling_distributions/seurat_resampling_at_scale.rds"))
resampling_res_processed <- process_undercover_result(resampling_res, sample_size_df)
saveRDS(object = resampling_res_processed, paste0(sceptre2_results_dir,
                                                  "resampling_distributions/seurat_resampling_at_scale_processed.rds"))

# pc results
rm(list = c("resampling_res", "resampling_res_processed"))
pc_res <- readRDS(paste0(sceptre2_results_dir, "positive_control_analysis/pc_results.rds"))
min_p <- pc_res |> filter(method == "sceptre", p_value > 0) |> pull(p_value) |> min()
pc_res_processed <- process_pc_result(pc_res, sample_size_df) |>
  mutate(p_value = ifelse(p_value <= 0, min_p, p_value))
saveRDS(object = pc_res_processed,
        file = paste0(sceptre2_results_dir, "positive_control_analysis/pc_results_processed.rds"))
