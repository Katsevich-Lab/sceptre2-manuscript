library(tidyverse)

sceptre2_results_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")

# load the sample size data frame
sample_size_df <- readRDS(paste0(sceptre2_results_dir, "dataset_sample_sizes/n_nonzero_cells_per_grna.rds"))
sample_size_df$dataset_concat |> unique()

eff_ss_df <- sample_size_df |>
  dplyr::filter(dataset_concat %in% c("frangieh/ifn_gamma/gene",
                                      "papalexi/eccite_screen/gene",
                                      "schraivogel/enhancer_screen_chr11/gene")) |>
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
  sample_n(size = n_ex_per_eff_ss)

# run SCEPTRE on these various pairs
