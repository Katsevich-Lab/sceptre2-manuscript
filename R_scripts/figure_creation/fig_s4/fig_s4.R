library(ggplot2)
library(katlabutils)
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), 
                            "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)

# directory with results
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")

# results of undercover analysis
undercover_res <- readRDS(paste0(result_dir,
                                "undercover_grna_analysis/undercover_result_grp_1_0523_processed.rds")) |>
  dplyr::filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF) |>
  dplyr::filter(!(Method %in% c("NB regression (no covariates)", 
                           "NB regression (w/ covariates)", 
                           "SCEPTRE (no covariates)"))) |>
  dplyr::mutate(Method = forcats::fct_relevel(Method, "SCEPTRE", after = Inf)) |>
  dplyr::filter(dataset %in% c("frangieh_ifn_gamma_gene", "papalexi_eccite_screen_gene",
                               "frangieh_control_gene", "simulated_experiment_1_gene"))

# for two of the datasets, cut the effective sample size into four intervals; then, plot each method, faceting by effective sample size
undercover_res_w_bin <- undercover_res |>
  dplyr::group_by(dataset) |>
  dplyr::mutate(n_nonzero_trt_bin = bin(n_nonzero_treatment, 4L))
my_methods <- c("KS test", "MAST", "MIMOSCA", "t-test", "Seurat-Wilcox", "Seurat-NB", "SCEPTRE")
my_values <- my_cols[names(my_cols) %in% my_methods]
my_y_breaks <- c(1e-0, 1e-3, 1e-6)
my_x_breaks <- c(1e-0, 1e-2, 1e-4)

p_frangieh_gamma <- ggplot(data = undercover_res_w_bin |>
                       dplyr::filter(dataset == "frangieh_ifn_gamma_gene"),
                     mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  stat_qq_band() +
  facet_grid(n_nonzero_trt_bin ~ dataset_rename) +
  scale_x_continuous(trans = revlog_trans(10), breaks = my_x_breaks) +
  scale_y_continuous(trans = revlog_trans(10), breaks = my_y_breaks) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") + 
  scale_color_manual(values = my_values) +
  my_theme_no_legend


p_frangieh_control <- ggplot(data = undercover_res_w_bin |>
                             dplyr::filter(dataset == "frangieh_control_gene"),
                           mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  stat_qq_band() +
  facet_grid(n_nonzero_trt_bin ~ dataset_rename) +
  scale_x_continuous(trans = revlog_trans(10), breaks = my_x_breaks) +
  scale_y_continuous(trans = revlog_trans(10), breaks = my_y_breaks) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") + 
  scale_color_manual(values = my_values) +
  my_theme_no_legend

p_simulated <- ggplot(data = undercover_res_w_bin |>
                        dplyr::filter(dataset == "simulated_experiment_1_gene"),
                      mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  stat_qq_band() +
  facet_grid(n_nonzero_trt_bin ~ dataset_rename) +
  scale_x_continuous(trans = revlog_trans(10), breaks = my_x_breaks) +
  scale_y_continuous(trans = revlog_trans(10), breaks = my_y_breaks) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") + 
  scale_color_manual(values = my_values) +
  my_theme_no_legend

p_papalexi <- ggplot(data = undercover_res_w_bin |>
                       dplyr::filter(dataset == "papalexi_eccite_screen_gene"),
                     mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  stat_qq_band() +
  facet_grid(n_nonzero_trt_bin ~ dataset_rename) +
  scale_x_continuous(trans = revlog_trans(10), breaks = my_x_breaks) +
  scale_y_continuous(trans = revlog_trans(10), breaks = my_y_breaks) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") + 
  scale_color_manual(values = my_values) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.margin = margin(t = 0, unit = "cm"),
        legend.key=element_blank()) +
  guides(color = guide_legend(
    override.aes = list(size = 2.5, fill = NA)))
legend <- cowplot::get_legend(p_papalexi)
p_papalexi <- p_papalexi + my_theme_no_legend + ylab("")
# create figures 1 and 2
p_combined_a <- cowplot::plot_grid(cowplot::plot_grid(p_frangieh_gamma, p_papalexi, ncol = 2),
                                   legend, rel_heights = c(0.93, 0.07), nrow = 2)
p_combined_b <- cowplot::plot_grid(cowplot::plot_grid(p_frangieh_control, p_simulated, ncol = 2),
                                   legend, rel_heights = c(0.93, 0.07), nrow = 2)
# save figures 1 and 2
to_save_fp_a <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                       "sceptre2-manuscript/R_scripts/figure_creation/fig_s4/fig_s4_a.png")
to_save_fp_b <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                       "sceptre2-manuscript/R_scripts/figure_creation/fig_s4/fig_s4_b.png")

ggsave(filename = to_save_fp_a, plot = p_combined_a, device = "png", scale = 1.1, width = 6.5, height = 7.0, dpi = 330)
ggsave(filename = to_save_fp_b, plot = p_combined_b, device = "png", scale = 1.1, width = 6.5, height = 7.0, dpi = 330)

