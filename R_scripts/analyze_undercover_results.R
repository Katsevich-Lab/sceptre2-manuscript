undercover_res_fps <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                             "results/undercover_grna_analysis/undercover_result_grp_size_", seq(1, 3), ".rds")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/figures/undercover_figs/")

# source the functions
funct_script <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                       "sceptre2-manuscript/R_scripts/analyze_undercover_results_plot_functs.R")
source(funct_script)

# load packages
library(ggplot2)
library(katlabutils)
library(future.apply)
plan(multisession)


future_lapply(X = seq(1, 3), FUN = function(i) {
  undercover_res <- readRDS(undercover_res_fps[i]) |> dplyr::mutate(clock_time = NULL, max_ram = NULL)
  undercover_res <- combine_schraivogel_enhancer_screens(undercover_res)
  if (!perform_sanity_check(undercover_res)) {
    stop("Sanity check failed.")
  }
  undercover_res <- update_dataset_names(undercover_res, TRUE)
  undercover_res_no_chrom <- undercover_res |>
    dplyr::filter(!(dataset %in% c("liscovitch_experiment_big_chromatin",
                                   "liscovitch_experiment_small_chromatin")))
  undercover_res_chrom_only <- undercover_res |>
    dplyr::filter(dataset %in% c("liscovitch_experiment_big_chromatin",
                                 "liscovitch_experiment_small_chromatin"))
  
  # 1. create the qq-plots
  trans_qq_plot <- make_trans_qq_plot(undercover_res_no_chrom)
  untrans_qq_plot <- make_untrans_qq_plot(undercover_res_no_chrom)
  lisc_plot <- make_trans_qq_plot(undercover_res_chrom_only)
  
  #ggsave(filename = paste0(fig_dir, "trans_qq_no_chrom_grp=", i, ".png"),
  #       plot = trans_qq_plot, device = "png", scale = 1.25, dpi = 330, width = 9, height = 6)
  #ggsave(filename = paste0(fig_dir, "untrans_qq_no_chrom_grp=", i, ".png"),
  #       plot = untrans_qq_plot, device = "png", scale = 1.25, dpi = 330, width = 9, height = 6)
  #ggsave(file = paste0(fig_dir, "trans_qq_chrom_only_grp=", i, ".png"),
  #       plot = lisc_plot, device = "png", scale = 1.4, dpi = 330, width = 4, height = 2)
  
  # 2. create n pairs rejected plots
  n_rejected_df <- compute_n_bonf_rejected(undercover_res = undercover_res_no_chrom, alpha = 0.05)
  n_rejected_plot <- make_n_rejected_pairs_plot(n_rejected_df = n_rejected_df, y_max = 1e5)
  ggsave(filename = paste0(fig_dir, "barplot_no_chrom_grp=", i, ".png"),
         plot = n_rejected_plot, device = "png", scale = 1.25, dpi = 330, width = 9, height = 6)
  
  # 3. create plot exploring connection between effective sample size and rejection rate
  
})
