library(ggplot2)

shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                            "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
correlated_sim_res <- readRDS(paste0(result_dir, "extra_analyses/correlated_sim_result.rds"))

z_star <- correlated_sim_res$resamp_dist$camp_star
p_val <- correlated_sim_res$out_m[1,"p_camp"]
density_df <- data.frame(z = correlated_sim_res$resamp_dist$camp_null)
p <- ggplot() + geom_histogram(aes(x = z, y = after_stat(density)),
                          data = density_df,  color = "black",
                          fill = "grey85", bins = 20) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = c(0.0, .01))) +
  theme(axis.ticks.x =  element_blank(),
        axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y =  element_blank()) +
  geom_vline(xintercept = z_star, col = "purple", linewidth = 1) +
  annotate("rect", xmin = z_star, xmax = Inf,
           ymin = 0, ymax = Inf, fill = "purple", alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax = -z_star,
           ymin = 0, ymax = Inf, fill = "purple", alpha = 0.15) +
  xlab("Null z-scores") +
  xlim(c(-3.5, 3.5))
ggsave(filename = paste0(.get_config_path("LOCAL_CODE_DIR"),
                         "sceptre2-manuscript/R_scripts/figure_creation/fig_3/histogram.png"),
       plot = p, device = "png", width = 4, height = 3)

# cowplot option

left <- cowplot::plot_grid(NULL, p + ggtitle("\n"), NULL, NULL)
right <- cowplot::plot_grid(NULL, NULL, ncol = 1, labels = c("b", "c"))
fig <- cowplot::plot_grid(left, right, ncol = 2, rel_widths = c(2/3, 1/3), labels = c("a"))  
ggsave(filename = paste0(.get_config_path("LOCAL_CODE_DIR"),
                         "sceptre2-manuscript/R_scripts/figure_creation/fig_3/cowplot.pdf"),
       plot = fig, device = "pdf", scale = 1.2, width = 7.5, height = 5)
