library(tidyverse)
library(katlabutils)
library(cowplot)

# Load scripts and results
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
extra_analyses_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/extra_analyses/")

grna_res <- readRDS(paste0(extra_analyses_dir, "papalexi_grna_confounding_tests.rds"))
gene_res <- readRDS(paste0(extra_analyses_dir, "papalex_gene_confounding_tests.rds"))

p1 <- ggplot(data = data.frame(p_value = grna_res),
             mapping = aes(y = p_value)) +
  stat_qq_points(ymin = 1e-8) +
  stat_qq_band() +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  my_theme +
  ggtitle("gRNA presence/absence")

p2 <- ggplot(data = data.frame(p_value = gene_res),
             mapping = aes(y = p_value)) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  my_theme +
  ggtitle("Gene expression")

p_out <- cowplot::plot_grid(p1, p2, nrow = 1, ncol = 2, labels = "auto")

to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_s5/fig_s5.png")
ggsave(filename = to_save_fp, plot = p_out, device = "png", scale = 1.2, width = 6, height = 2.5, dpi = 330)
