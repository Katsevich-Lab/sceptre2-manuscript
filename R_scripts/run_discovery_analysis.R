#!/usr/bin/env Rscript
# The purpose of this script is to: (1) run the (approximate) calibration check for Papalexi and Frangieh IFN-gamma and (2) run the trans discovery analysis for these two datasets. We time and profile the memory of all four analyses from within Linux.

args <- commandArgs(trailingOnly = TRUE)

dataset <- args[1] # "papalexi" or "frangieh" or ""
analysis_type <- args[2] # "calibration" or "discovery"

cat(dataset); cat("\n")
cat(analysis_type); cat("\n")

library(sceptre)
library(Matrix)

LOCAL_SCEPTRE2_DATA_DIR <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
if (dataset == "papalexi") {
  objects_fp <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "data/papalexi/eccite_screen/r_objects.rds")
  f_name <- paste0("papalexi_gene_", analysis_type, "_res.rds")
} else {
  objects_fp <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "data/frangieh/ifn_gamma/r_objects.rds")
  f_name <- paste0("frangieh_ifn_gamma_", analysis_type, "_res.rds")
}

calibration_check <- analysis_type == "calibration"
l <- readRDS(objects_fp)
gc() |> invisible()

res <- run_sceptre_lowmoi(response_matrix = l$response_matrix,
                          grna_matrix = l$grna_matrix,
                          covariate_data_frame = l$covariate_data_frame,
                          grna_group_data_frame = l$grna_group_data_frame,
                          formula_object = l$formula_object,
                          response_grna_group_pairs = l$response_grna_group_pairs,
                          calibration_check = calibration_check)

gc() |> invisible()


result_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "results/discovery_analyses/")
saveRDS(object = res, file = paste0(result_dir, f_name))