library(tidyverse)
sceptre2_results_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")

# load the sample size data frame
sample_size_df <- readRDS(paste0(sceptre2_results_dir, "dataset_sample_sizes/n_nonzero_cells_per_grna.rds"))
sample_size_df$dataset_concat |> unique()
my_datasets <- c("frangieh/ifn_gamma/gene",
                 "papalexi/eccite_screen/gene",
                 "schraivogel/enhancer_screen_chr11/gene")

eff_ss_df <- sample_size_df |>
  dplyr::filter(dataset_concat %in% my_datasets) |>
  group_by(feature_id, dataset_concat, target) |>
  summarize(eff_ss = sum(n_nonzero_cells)) |>
  ungroup()

# set different effective sample size levels; sample from each
eff_ss_levels <- c(7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30)
n_ex_per_eff_ss <- 10
set.seed(4)
my_samp <- eff_ss_df |> 
  filter(eff_ss %in% eff_ss_levels) |>
  group_by(dataset_concat, eff_ss) |>
  sample_n(size = n_ex_per_eff_ss) |>
  rename(response_id = feature_id, grna_group = target)
rm(sample_size_df, eff_ss_df); gc()

to_save <- lapply(X = my_datasets, FUN = function(my_dataset) {
  response_grna_group_pairs <- my_samp |> filter(dataset_concat == !!my_dataset)
  response_odm <- lowmoi::load_dataset_modality(my_dataset)
  grna_dataset_name <- lowmoi::get_grna_dataset_name(my_dataset, "assignment")
  grna_odm <- lowmoi::load_dataset_modality(grna_dataset_name)
  mm_odm <- ondisc::multimodal_ondisc_matrix(list(response = response_odm, grna = grna_odm)) |>
    lowmoi::process_multimodal_odm()
  res <- sceptre2::run_sceptre_low_moi(mm_odm = mm_odm,
                                       response_grna_group_pairs = response_grna_group_pairs,
                                       form = response_odm@misc$sceptre_formula,
                                       response_modality_name = "response",
                                       grna_modality_name = "grna",
                                       grna_group_column_name = "target",
                                       B = 500000,
                                       max_b_per_batch = 250000,
                                       in_memory = TRUE,
                                       statistic = "full",
                                       return_dist = TRUE)
  
  res_order <- left_join(res |> select(response_id, grna_group),
                         response_grna_group_pairs,
                         by = c("response_id", "grna_group"))
  res$eff_ss <- res_order$eff_ss
  res <- res |>
    mutate(dataset = factor(my_dataset)) |>
    relocate(response_id, grna_group, eff_ss, dataset)
  res
}) |> data.table::rbindlist()

# save
saveRDS(object = to_save, file = paste0(sceptre2_results_dir, "resampling_distributions/sceptre_resampling_dists.rds"))

# AN EXPLORATION OF USING 
if (FALSE) {
  library(evmix)
  # load empirical distributions
  to_save <- readRDS(paste0(sceptre2_results_dir, "resampling_distributions/sceptre_resampling_dists.rds"))
  
  # index of to_save to investigate
  i <- sample(seq(1, nrow(to_save)), 1)
  to_save[i, 1:4] # response, grna, effective sample size (i.e., N nonzero treatment cells, dataset)
  
  # investigate left tail or right tail
  left_tail <- FALSE
  z_null <- (if (left_tail) -1 else 1) * (to_save[i, 5:ncol(to_save_sub)] |> as.numeric())
  
  # SELECT CUTTOF; EXAMINE TAIL DIST
  u <- 2.5
  z_null_tail <- z_null[z_null > u]
  max(z_null_tail)
  # histogram
  hist(z_null_tail, freq = FALSE)
  # fit GDP via MLE
  fit <- fgpd(x = z_null, u = u)
  # plot fitted dist
  x_grid <- seq(u, 6, length.out = 1000)
  y_dens <- dgpd(x = xgrid, u = 2, sigmau = fit$mle[1], xi = fit$mle[2])
  lines(x = x_grid, y = y_dens, pch = 19, cex = 0.4, col = "red")
  # compare empirical (conditional) p-value to fitted p-value
  z_star <- 3.75
  emp_p <- mean(z_null > u) *  (1 - sceptre2:::compute_empirical_p_value(z_star, z_null_tail, side = "left"))
  p_p <- mean(z_null > u) * (1 - pgpd(q = z_star, u = u, sigmau = fit$mle[1], xi = fit$mle[2]))
  emp_p
  p_p
  emp_p/p_p
}
