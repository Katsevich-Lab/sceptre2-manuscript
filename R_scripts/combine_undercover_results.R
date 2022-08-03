# combines undercover results undercover_result_grp_size_n.rds and undercover_result_grp_size_n_new_perm.rds, replacing
# permutation_test results of the former with those of the latter

n <- 1
undercover_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/undercover_grna_analysis/")
res <- paste0(undercover_dir, "undercover_result_grp_size_", n ,".rds") |> readRDS()
res_new_perm <- paste0(undercover_dir, "undercover_result_grp_size_", n, "_new_perm.rds") |> readRDS()
res_new <- rbind(res |> dplyr::filter(method != "permutation_test"), res_new_perm)
saveRDS(object = res_new, file = paste0(undercover_dir, "undercover_result_grp_size_", n ,".rds"))

# Combines positive control results
pc_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/positive_control_analysis/")
r1 <- readRDS(paste0(pc_dir, "pc_result.rds"))
r2 <- readRDS(paste0(pc_dir, "pc_results_extra.rds"))

out <- rbind(dplyr::filter(r1, dataset != "schraivogel/ground_truth_tapseq/gene"),
             dplyr::distinct(r2))

saveRDS(object = out, file = paste0(pc_dir, "pc_result.rds"))
