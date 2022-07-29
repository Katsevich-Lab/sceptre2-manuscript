undercover_res_fps <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                             "results/undercover_grna_analysis/undercover_result_grp_size_", seq(1, 3), ".rds")
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/figures/undercover_figs/")

# source the functions
funct_script <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                       "sceptre2-manuscript/writeups/digging_into_undercover/analyze_undercover_results_plot_functs.R")

# load packages and script
source(funct_script)
library(katlabutils)
library(tidyverse)
#library(future.apply)
#plan(multisession)


lapply(X = seq(1, 3), FUN = function(i) {
  undercover_res <- readRDS(undercover_res_fps[i]) |>
    dplyr::mutate(clock_time = NULL, max_ram = NULL)|>
    dplyr::mutate(dataset_slash = replace_dataset_underscore_with_slash(dataset))
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
  if (FALSE) {
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
  }
  
  # 2. create n pairs rejected plots
  n_rejected_df <- compute_n_bonf_rejected(undercover_res = undercover_res_no_chrom, alpha = 0.05)
  n_rejected_plot <- make_n_rejected_pairs_plot(n_rejected_df = n_rejected_df, y_max = 1e5)
  
  if (FALSE) {
  ggsave(filename = paste0(fig_dir, "barplot_no_chrom_grp=", i, ".png"),
         plot = n_rejected_plot, device = "png", scale = 1.25, dpi = 330, width = 9, height = 6)
  }

  # 3. create histograms for the permutation, seurat de, and NB regression methods
  methods <- c("permutation_test", "seurat_de", "nb_regression")
  if (FALSE) {
  for (j in seq(1, length(methods))) {
    method <- methods[j]
    fill <- my_cols[j]
    hist_p <- make_histograms(undercover_res |> filter(method == !!method),
                              fill = my_cols[j])
    ggsave(filename = paste0(fig_dir, "histogram_", method, "_grp=", i, ".png"),
           plot = hist_p, device = "png", scale = 1.1, width = 8, height = 5, dpi = 330)
  }
  }
  
  # 4. create plot exploring connection between effective sample size and rejection rate
  undercover_res <- append_n_undercover_cells(undercover_res)
  n_rejected_n_cells_df <- get_n_reject_n_undercover_cells_df(undercover_res, alpha = 0.1)
  
  methods <- c("seurat_de", "nb_regression")
  for (j in seq(1, length(methods))) {
    method <- methods[j]
    effective_samp_size_n_rejected_assoc_plot <- associate_n_undercover_cells_w_n_rejections(n_rejected_n_cells_df |> filter(method == !!method))
    ggsave(filename = paste0(fig_dir, "assoc_plot_", method, "_grp=", i, ".png"),
           plot = effective_samp_size_n_rejected_assoc_plot,
           device = "png", scale = 1.1, width = 8, height = 5,
           dpi = 330)
  }
})


# investigate outlier: non-targeting-00022 on schraivogel_ground_truth_perturbseq_gene
res_grp_1 <- readRDS(undercover_res_fps[1]) |>
  filter(dataset == "schraivogel_ground_truth_perturbseq_gene")
p1 <- res_grp_1 |>
  update_dataset_names() |>
  make_trans_qq_plot()
p2 <- res_grp_1 |>
  filter(undercover_grna != "non-targeting-00022") |>
  update_dataset_names() |>
  make_trans_qq_plot()
ggsave(filename = paste0(fig_dir, "schraiv_perturb_qq_trans_outlier_removed.png"),
       plot = p1, device = "png", scale = 1.1, width = 5, height = 5, dpi = 330)
ggsave(filename = paste0(fig_dir, "schraiv_perturb_qq_trans_outlier_retained.png"),
       plot = p2, device = "png", scale = 1.1, width = 5, height = 5, dpi = 330)


# make NTC-specific qq-plots
# first, compute the number of unique NTCs in each dataset
undercover_res <- readRDS(undercover_res_fps[1]) |>
  dplyr::mutate(clock_time = NULL, max_ram = NULL)|>
  dplyr::mutate(dataset_slash = replace_dataset_underscore_with_slash(dataset)) |>
  combine_schraivogel_enhancer_screens() |>
  update_dataset_names(TRUE)

# loop over select methods and datasets, creating the gRNA-specific qq-plot for each
methods <- c("seurat_de", "nb_regression", "permutation_test")
datasets <- as.character(unique(undercover_res$dataset))
datasets <- datasets[!(datasets %in% c("papalexi_eccite_screen_protein",
                                       "liscovitch_experiment_big_chromatin",
                                       "liscovitch_experiment_small_chromatin"))]

for (curr_method in methods) {
  for (curr_dataset in datasets) {
    print(paste0("Working on ", curr_dataset, " and ", curr_method))
    undercover_res_sub <- undercover_res |>
      filter(method == curr_method & dataset == curr_dataset)
    if (curr_dataset %in% c("frangieh_co_culture_gene",
                            "frangieh_control_gene",
                            "frangieh_ifn_gamma_gene")) {
      undercover_res_sub <- undercover_res_sub |> filter(undercover_grna != "NO-SITE-707")
      scale <- 1.5
    }
    p <- make_grna_specific_plots(undercover_res_sub)
    # save
    ggsave(filename = paste0(fig_dir, "grna_wise_qqplots_", curr_dataset, "_", curr_method, ".png"),
           plot = p, device = "png", scale = 1.5, width = 14, height = 8, dpi = 330)
  }
}
