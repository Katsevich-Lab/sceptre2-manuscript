library(lowmoi)
library(tidyverse)
conflicts_prefer(dplyr::filter)
sceptre2_results_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
sample_size_df <- readRDS(paste0(sceptre2_results_dir, "dataset_sample_sizes/n_nonzero_cells_per_grna.rds"))

# undercover res grp = 1
undercover_res <- readRDS(paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1_0423.rds"))
sceptre_undercover_res <- readRDS(paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1_sceptre_0423.rds"))
undercover_res <- rbind(undercover_res, sceptre_undercover_res)
undercover_res_processed <- process_undercover_result(undercover_res, sample_size_df)
saveRDS(object = undercover_res_processed, 
        paste0(sceptre2_results_dir, "undercover_grna_analysis/undercover_result_grp_1_0423_processed.rds"))

# resampling results
resampling_res <- readRDS(paste0(sceptre2_results_dir,
                                 "resampling_distributions/seurat_resampling_at_scale.rds"))
resampling_res_processed <- process_undercover_result(resampling_res, sample_size_df)
saveRDS(object = resampling_res_processed, paste0(sceptre2_results_dir,
                                                  "resampling_distributions/seurat_resampling_at_scale_processed.rds"))

# pc result
pc_res <- readRDS(paste0(sceptre2_results_dir, "positive_control_analysis/pc_results_0423.rds"))
sceptre_pc_res <- readRDS(paste0(sceptre2_results_dir, "positive_control_analysis/pc_results_sceptre_0423.rds"))
pc_res <- rbind(pc_res, sceptre_pc_res)
pc_res_processed <- process_pc_result(pc_res, sample_size_df)
saveRDS(object = pc_res_processed,
        file = paste0(sceptre2_results_dir, "positive_control_analysis/pc_results_processed.rds"))

# discovery results
disc_res <- readRDS(paste0(sceptre2_results_dir, "discovery_analyses/discovery_results_0423.rds"))
sceptre_disc_res <- readRDS(paste0(sceptre2_results_dir, "discovery_analyses/disc_analysis_sceptre_0423.rds"))
disc_res <- rbind(disc_res, sceptre_disc_res)
disc_res_processed <- process_pc_result(pc_res = disc_res,
                                        sample_size_df = sample_size_df)
saveRDS(object = disc_res_processed,
        file = paste0(sceptre2_results_dir, "discovery_analyses/discovery_results_0423_processed.rds"))
