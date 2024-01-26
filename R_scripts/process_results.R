library(lowmoi)
library(tidyverse)
conflicts_prefer(dplyr::filter)
sceptre2_results_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
sample_size_df <- readRDS(paste0(sceptre2_results_dir, "dataset_sample_sizes/n_nonzero_cells_per_grna.rds"))

# undercover res grp = 1
undercover_res <- readRDS(paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1_0523.rds"))
undercover_res_processed <- process_undercover_result(undercover_res, sample_size_df)
saveRDS(object = undercover_res_processed, 
        paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1_0523_processed.rds"))

undercover_res_extra <- readRDS(paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1_extras_0523.rds"))
undercover_res_extra_processed <- suppressWarnings(process_undercover_result(undercover_res_extra, sample_size_df))
saveRDS(object = undercover_res_extra_processed,
        paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1_extras_0523_processed.rds"))

# pc result
pc_res <- readRDS(paste0(sceptre2_results_dir, "positive_control_analysis/pc_results_0124.rds"))
pc_res_processed <- suppressWarnings(process_pc_result(pc_res, sample_size_df))
saveRDS(object = pc_res_processed,
        file = paste0(sceptre2_results_dir, "positive_control_analysis/pc_results_0124_processed.rds"))

pc_res_sceptre_unfiltered <- readRDS(paste0(sceptre2_results_dir, "positive_control_analysis/pc_results_sceptre_unfiltered_0523.rds"))
pc_res_sceptre_unfiltered_processed <- suppressWarnings(process_pc_result(pc_res_sceptre_unfiltered, sample_size_df))
saveRDS(object = pc_res_sceptre_unfiltered_processed,
        file = paste0(sceptre2_results_dir, "positive_control_analysis/pc_results_sceptre_unfiltered_0523_processed.rds"))


# resampling results
resampling_res <- readRDS(paste0(sceptre2_results_dir,
                                 "resampling_distributions/seurat_resampling_at_scale.rds"))
resampling_res_processed <- suppressWarnings(process_undercover_result(resampling_res, sample_size_df))
saveRDS(object = resampling_res_processed, paste0(sceptre2_results_dir,
                                                  "resampling_distributions/seurat_resampling_at_scale_processed.rds"))

# discovery result
discovery_res <- readRDS(paste0(sceptre2_results_dir, "discovery_analyses/discovery_results_0423.rds"))
discovery_res_processed <- suppressWarnings(process_pc_result(discovery_res, sample_size_df))
saveRDS(discovery_res_processed, paste0(sceptre2_results_dir, "discovery_analyses/discovery_results_0423_processed.rds"))
