# Load packages
library(tidyverse)
library(katlabutils)
library(cowplot)
conflicts_prefer(dplyr::filter)

# Load scripts and results
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
undercover_res_with_biorep <- readRDS(paste0(result_dir, "undercover_grna_analysis/undercover_result_grp_1_0523_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF, n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF) |>
  filter(Method == "KS test", dataset == "papalexi_eccite_screen_gene") |>
  mutate(bio_rep = "With bio. rep")
undercover_res_without_biorep <- readRDS(paste0(result_dir, "undercover_grna_analysis/undercover_result_grp_1_extras_0523_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF, n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF) |>
  filter(Method == "KS test", dataset == "papalexi_eccite_screen_gene") |>
  mutate(bio_rep = "Without bio. rep")
res <- rbind(undercover_res_with_biorep, undercover_res_without_biorep)

p1 <- ggplot(data = res, mapping = aes(y = p_value, col = bio_rep)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  stat_qq_band() +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") + my_theme_no_legend +
  scale_color_manual(values = c("purple", "purple4"))

p2 <- ggplot(data = res, mapping = aes(y = p_value, col = bio_rep)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4)) +
  scale_y_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4)) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") + my_theme +
  scale_color_manual(values = c("purple", "purple4")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.75, 0.1),
        legend.text=element_text(size = 8.0),
        legend.margin=margin(t = 0, unit='cm'),
        legend.background = element_blank()) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.1,
    default.unit="inch",
    override.aes = list(size = 2.5)))

p <- plot_grid(p1, p2, labels = "auto")
to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/fig_s4/fig_s4.png")
ggsave(filename = to_save_fp, plot = p, device = "png", scale = 1, width = 6, height = 2.5, dpi = 330)
