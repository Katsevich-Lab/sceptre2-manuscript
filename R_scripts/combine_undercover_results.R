# combines undercover results undercover_result_grp_size_n.rds and undercover_result_grp_size_n_new_perm.rds, replacing
# permutation_test results of the former with those of the latter

n <- 1
undercover_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/undercover_grna_analysis/")
res <- paste0(undercover_dir, "undercover_result_grp_size_", n ,".rds") |> readRDS()
res_new_perm <- paste0(undercover_dir, "undercover_result_grp_size_", n, "_new_perm.rds") |> readRDS()
res_new <- rbind(res |> dplyr::filter(method != "permutation_test"), res_new_perm)
saveRDS(object = res_new, file = paste0(undercover_dir, "undercover_result_grp_size_", n ,".rds"))
