library(tidyverse)
library(katlabutils)
library(cowplot)

# load functions and data
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/extra_analyses/")
correlated_res <- readRDS(paste0(result_dir, "correlated_sim_result.rds"))
uncorrelated_res <- readRDS(paste0(result_dir, "uncorrelated_sim_result.rds"))

p1 <- as.data.frame(correlated_res$out_m) |>
  tidyr::pivot_longer(cols = c("p_theory", "p_camp", "p_perm"),
                      names_to = "Method", values_to = "p_value") |>
  mutate(Method = fct_recode(Method, "NB Regression" = "p_theory",
                             "SCEPTRE" = "p_camp",
                             "Permutation test" = "p_perm")) |>
  ggplot(mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  my_theme +
  theme(legend.position = c(0.75, 0.17),
        legend.title = element_blank()) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.1,
    default.unit="inch")) +
  ggtitle("Treatment confounded,\nGLM specified correctly") +
  scale_color_manual(values = my_cols)


p2 <- as.data.frame(uncorrelated_res$out_m) |>
  tidyr::pivot_longer(cols = c("p_theory", "p_camp", "p_perm"),
                      names_to = "Method", values_to = "p_value") |>
  mutate(Method = fct_recode(Method, "NB Regression" = "p_theory",
                             "SCEPTRE" = "p_camp",
                             "Permutation test" = "p_perm")) |>
  ggplot(mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  my_theme +
  theme(legend.position = "none") +
  ggtitle("Treatment unconfounded,\nGLM specified incorrectly") +
  scale_color_manual(values = my_cols)

z_grid <- seq(-4, 4, length.out = 1000)
density_df <- data.frame(density = dnorm(z_grid),
                         z_grid = z_grid)
histogram_df <- data.frame(z_null = correlated_res$resamp_dist$camp_null)
p3 <- ggplot() + geom_histogram(aes(x = z_null, y = after_stat(density)),
                          data = histogram_df,
                          boundary = 0,
                          fill = "grey85",
                          color = "black",
                          bins = 25) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = c(0.0, .01))) +
  geom_line(aes(x = z_grid, y = density, col = "N(0,1) density"),
            data = density_df, linewidth = 0.7) +
  scale_color_manual(values = c("N(0,1) density" = "purple")) +
  xlab("z null") +
  ylab("") +
  theme(legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(),
        legend.margin=margin(t = -0.5, unit='cm')) +
  ggtitle("SCEPTRE null z-scores") +
  xlim(-4, 4)
  

histogram_df <- data.frame(z_null = uncorrelated_res$resamp_dist$camp_null)
p4 <- ggplot() + geom_histogram(aes(x = z_null, y = after_stat(density)),
                                data = histogram_df,
                                boundary = 0,
                                fill = "grey85",
                                color = "black",
                                bins = 25) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = c(0.0, .01))) +
  geom_line(aes(x = z_grid, y = density),
            data = density_df, linewidth = 0.7, col = "purple") +
  xlab("z null") +
  ylab("") +
  ggtitle("SCEPTRE null z-scores") +
  xlim(-6, 6)

fig <- plot_grid(p1, p2, p3, p4, nrow = 2, rel_heights = c(0.55, 0.45),
                 labels = c("a", "c", "b", "d"), byrow = TRUE)

to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_s4/fig_s4.pdf")
ggsave(filename = to_save_fp, plot = fig, device = "pdf", scale = 1, width = 6.5, height = 5)
