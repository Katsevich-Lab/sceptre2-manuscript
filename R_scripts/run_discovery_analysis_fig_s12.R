#!/usr/bin/env Rscript
# The purpose of this script is to: (1) run the calibration check for Papalexi and Frangieh IFN-gamma and (2) run the trans discovery analysis for these two datasets. We time and profile the memory of all four analyses from within Linux.

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1] # "papalexi" or "frangieh"
full_statistic <- as.logical(args[2]) # TRUE (for full) or FALSE (for residuals)

cat(paste0("dataset: ", dataset, "\n"))
cat(paste0("full statistic: ", full_statistic, "\n"))

library(sceptre)
library(Matrix)

LOCAL_SCEPTRE2_DATA_DIR <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
if (dataset == "papalexi") {
  objects_fp <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "data/papalexi/eccite_screen/r_objects.rds")
} else {
  objects_fp <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "data/frangieh/control/r_objects.rds")
}
l <- readRDS(objects_fp)
gc() |> invisible()

###################################################
# Prepare the analysis by creating a sceptre object
###################################################
# import data
response_matrix <- l$response_matrix
grna_matrix <- l$grna_matrix
grna_target_data_frame <- l$grna_group_data_frame
sceptre_object <- import_data(response_matrix = l$response_matrix,
                              grna_matrix = l$grna_matrix,
                              grna_target_data_frame = l$grna_group_data_frame |>
                                dplyr::rename(grna_target = grna_group),
                              moi = "low")

##################
# Run the analysis
##################
trans_pairs <- construct_trans_pairs(sceptre_object)
pc_pairs <- construct_positive_control_pairs(sceptre_object)
sceptre_object <- set_analysis_parameters(sceptre_object,
                                          discovery_pairs = trans_pairs,
                                          positive_control_pairs = pc_pairs,
                                          full_test_stat = full_statistic)
# assign grnas
sceptre_object <- sceptre_object |> assign_grnas()
sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = 0.1)

# run the calibration check
calibration_check_time <- system.time(
  sceptre_object <- sceptre_object |>
    run_calibration_check(parallel = TRUE, n_processors = 8) 
)[["elapsed"]]
# run the positive control analysis
sceptre_object <- sceptre_object |> run_power_check(parallel = TRUE, n_processors = 8)
# run the discovery analysis
discovery_analysis_time <- system.time(
  sceptre_object <- sceptre_object |>
    run_discovery_analysis(parallel = TRUE, n_processors = 8)
)[["elapsed"]]

# obtain the calibration check and discovery analysis results
dir_name <- paste0(dataset, "_", (if (full_statistic) "full_stat" else "resid_stat"))
write_outputs_to_directory(sceptre_object,
                           directory = paste0(LOCAL_SCEPTRE2_DATA_DIR, "results/discovery_analyses/fig_s12/", dir_name))
running_times <- c(calibration_check_time = calibration_check_time,
  discovery_analysis_time = discovery_analysis_time)
saveRDS(running_times, file = paste0(LOCAL_SCEPTRE2_DATA_DIR, "results/discovery_analyses/fig_s12/", dir_name, "/running_times.rds"))
