library(lowmoi)
sceptre2_results_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
undercover_res <- readRDS(paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1.rds"))
undercover_res_extra <- readRDS(paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1_debug.rds"))
undercover_res <- rbind(undercover_res, undercover_res_extra)
sample_size_df <- readRDS(paste0(sceptre2_results_dir, "dataset_sample_sizes/n_nonzero_cells_per_grna.rds"))
undercover_res_processed <- process_undercover_result(undercover_res, sample_size_df) |>
  mutate(p_value = ifelse(p_value <= 0, 1e-8, p_value))
saveRDS(undercover_res_processed, paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1_processed.rds"))
