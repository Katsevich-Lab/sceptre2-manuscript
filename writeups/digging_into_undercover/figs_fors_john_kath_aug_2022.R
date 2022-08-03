undercover_res_fps <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                             "results/undercover_grna_analysis/undercover_result_grp_size_", seq(1, 3), ".rds")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/writeups/digging_into_undercover/figs_for_john_kat_aug_2022/")
if (!dir.exists(fig_dir)) dir.create(fig_dir)

# source the functions
funct_script <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                       "sceptre2-manuscript/writeups/digging_into_undercover/analyze_undercover_results_plot_functs.R")

# load packages and script
source(funct_script)
library(katlabutils)
library(tidyverse)
library(future.apply)

plan(multisession)
future_lapply(X = seq(1, 3), FUN = function(i) {
  undercover_res <- readRDS(undercover_res_fps[i]) |>
    dplyr::mutate(clock_time = NULL, max_ram = NULL)|>
    dplyr::filter(method %in% c("seurat_de", "liscovitch_method", "schraivogel_method", "weissman_method", "mimosca")) |>
    combine_schraivogel_enhancer_screens() |>
    update_dataset_names(TRUE)
  
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
    
    ggsave(filename = paste0(fig_dir, "trans_qq_no_chrom_grp=", i, ".png"),
           plot = trans_qq_plot, device = "png", scale = 1.25, dpi = 330, width = 9, height = 6)
    ggsave(filename = paste0(fig_dir, "untrans_qq_no_chrom_grp=", i, ".png"),
           plot = untrans_qq_plot, device = "png", scale = 1.25, dpi = 330, width = 9, height = 6)
    ggsave(file = paste0(fig_dir, "trans_qq_chrom_only_grp=", i, ".png"),
           plot = lisc_plot, device = "png", scale = 1.4, dpi = 330, width = 4, height = 2)

  # 2. create n pairs rejected plots
  n_rejected_df <- compute_n_bonf_rejected(undercover_res = undercover_res_no_chrom, alpha = 0.05)
  n_rejected_plot <- make_n_rejected_pairs_plot(n_rejected_df = n_rejected_df, y_max = 1e5)
  
  ggsave(filename = paste0(fig_dir, "barplot_no_chrom_grp=", i, ".png"),
         plot = n_rejected_plot, device = "png", scale = 1.25, dpi = 330, width = 9, height = 6)
})

# Positive controls
# load data; set figure dir
pc_res_fp <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/positive_control_analysis/pc_result.rds")
pc_res <- readRDS(pc_res_fp) |>
  filter(method %in% c("seurat_de", "liscovitch_method", "schraivogel_method", "weissman_method", "mimosca")) |>
  arrange(dataset, method)
funct_script <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                       "sceptre2-manuscript/writeups/digging_into_undercover/analyze_undercover_results_plot_functs.R")
pc_res |>
  group_by(dataset, method) |>
  summarize(n_pairs = n()) |>
  summarize(same_n_pairs = all(diff(n_pairs) == 0)) |>
  pull(same_n_pairs) |>
  all()
pc_res <- pc_res |>
  dplyr::mutate(dataset = gsub(pattern = "/", replacement = "_", fixed = TRUE, x = dataset))
pc_res <- combine_schraivogel_enhancer_screens(pc_res) |>
  update_dataset_names(TRUE)
n_bonf_reject <- compute_n_bonf_rejected(pc_res)
p <- make_n_rejected_pairs_plot(n_rejected_df = n_bonf_reject, y_max = NULL, scales = "free", log_trans = FALSE)
ggsave(filename = paste0(fig_dir, "pc_results.png"), plot = p, device = "png", scale = 1.1, width = 8, height = 6, dpi = 330)
