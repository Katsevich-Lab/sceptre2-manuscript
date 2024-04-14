library(tidyverse)
library(katlabutils)
library(cowplot)

# load functions and data
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/extra_analyses/")
sim_res <- readRDS(paste0(result_dir, "simulation_study_res.rds"))

# p1: no problems (no confounding, correct model specification)
make_qq_plot <- function(result_matrix, tit, w_legend = FALSE) {
  p0 <- as.data.frame(result_matrix) |>
    tidyr::pivot_longer(cols = c("p_theory", "p_camp", "p_perm"),
                        names_to = "Method", values_to = "p_value") |>
    mutate(Method = fct_recode(Method, "NB regression (w/ covariates)" = "p_theory",
                               "SCEPTRE" = "p_camp",
                               "Permutation test" = "p_perm")) |>
    mutate(Method = fct_relevel(Method, "SCEPTRE", after = Inf)) |>
    ggplot(mapping = aes(y = p_value, col = Method)) +
    stat_qq_points(ymin = 1e-8, size = 0.9) +
    scale_x_reverse() +
    scale_y_reverse() +
    labs(x = "Expected null p-value", y = "Observed p-value") +
    geom_abline(col = "black") +
    my_theme +
    guides(color = guide_legend(
      keywidth = 0.0,
      keyheight = 0.2,
      default.unit="inch",
      override.aes = list(size = 2.5))) +
    ggtitle(tit) +
    scale_color_manual(values = my_cols)
  
  if (w_legend) {
    p <- p0 + theme(legend.position = c(0.7, 0.15),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.35, 'cm'),
                    legend.margin = margin(t = -0.5, unit = 'cm')
                    )
  } else {
    p <- p0 + theme(legend.position = "none")
  }
  return(p)
}

p_a <- make_qq_plot(result_matrix = sim_res$sim_res_uncorrelated_correct_model$out_m,
                    tit = "Treatment unconfounded,\nGLM specified correctly", w_legend = TRUE)

p_b <- make_qq_plot(result_matrix = sim_res$sim_res_uncorrelated_misspec$out_m,
                    tit = "Treatment unconfounded,\nGLM specified incorrectly", w_legend = FALSE)

p_c <- make_qq_plot(result_matrix = sim_res$sim_res_correlated_corret_model$out_m,
                    tit = "Treatment confounded,\nGLM specified correctly", w_legend = FALSE)

p_d <- make_qq_plot(result_matrix = sim_res$sim_res_correlated_misspec$out_m,
                    tit = "Treatment confounded,\nGLM specified incorrectly", w_legend = FALSE)

fig <- plot_grid(p_a, p_b, p_c, p_d, nrow = 2,
                 labels = "auto")

to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_s7/fig_s7.png")
ggsave(filename = to_save_fp, plot = fig, device = "png", scale = 1.1, width = 6.5, height = 5.5, dpi = 330)
