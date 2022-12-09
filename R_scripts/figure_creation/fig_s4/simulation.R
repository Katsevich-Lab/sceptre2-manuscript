library(tidyverse)
library(katlabutils)
library(cowplot)

# load functions and data
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/extra_analyses/")
correlated_res <- readRDS(paste0(result_dir, "correlated_sim_result.rds"))
uncorrelated_res <- readRDS(paste0(result_dir, "uncorrelated_sim_result.rds"))

#############################
# QQ-plot for confounded case
#############################
p1 <- as.data.frame(correlated_res$out_m) |>
  tidyr::pivot_longer(cols = c("p_theory", "p_camp", "p_perm"),
                      names_to = "Method", values_to = "p_value") |>
  mutate(Method = fct_recode(Method, "NB Regression" = "p_theory",
                             "SCEPTRE" = "p_camp",
                             "Vanilla permutation\ntest" = "p_perm")) |>
  ggplot(mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  my_theme +
  theme(legend.position = c(0.8, 0.2),
        legend.title = element_blank())


p2 <- as.data.frame(uncorrelated_res$out_m) |>
  tidyr::pivot_longer(cols = c("p_theory", "p_camp", "p_perm"),
                      names_to = "method", values_to = "p_value") |>
  ggplot(mapping = aes(y = p_value, col = method)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  my_theme +
  theme(legend.position = "none")


fig <- plot_grid(p1, p2, nrow = 1)
