# load data; set figure dir
pc_res_fp <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/positive_control_analysis/pc_result.rds")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/figures/undercover_figs/")

# load packages and data
library(tidyverse)
pc_res <- readRDS(pc_res_fp) |>
  arrange(dataset, method)

# load aux functions
funct_script <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                       "sceptre2-manuscript/writeups/digging_into_undercover/analyze_undercover_results_plot_functs.R")
source(funct_script)

# sanity check: verify that the number of p-values across datasets coincides
pc_res |>
  group_by(dataset, method) |>
  summarize(n_pairs = n()) |>
  summarize(same_n_pairs = all(diff(n_pairs) == 0)) |>
  pull(same_n_pairs) |>
  all()

# combine Schraivogel screens
pc_res <- pc_res |>
  dplyr::mutate(dataset = gsub(pattern = "/", replacement = "_", fixed = TRUE, x = dataset))
pc_res <- combine_schraivogel_enhancer_screens(pc_res) |>
  update_dataset_names(TRUE)

# compute number of Bonferoni rejections
n_bonf_reject <- compute_n_bonf_rejected(pc_res)

# plot the number of Bonferoni rejections
make_n_rejected_pairs_plot(n_rejected_df = n_bonf_reject, y_max = NULL, scales = "free", log_trans = FALSE)
