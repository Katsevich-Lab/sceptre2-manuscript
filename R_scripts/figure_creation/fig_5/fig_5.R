library(sceptre)
library(ggplot2)
fig_5_new_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/fig_5/")
papa_result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/discovery_analyses/")
calibration_result <- readRDS(paste0(papa_result_dir, "papalexi_gene_calibration_res.rds"))
discovery_result <- readRDS(paste0(papa_result_dir, "papalexi_gene_discovery_res.rds"))

# plot calibration result
calibration_plots <- plot_calibration_result(calibration_result = calibration_result, return_indiv_plots = TRUE)
transformed_qq_plot <- calibration_plots[[2]]
lfc_plot <- calibration_plots[[3]]
ggsave(plot = transformed_qq_plot, filename = paste0(fig_5_new_dir, "transformed_qq_plot.png"),
       device = "png", scale = 1.2, width = 2.5, height = 2, dpi = 330)
ggsave(plot = lfc_plot, filename = paste0(fig_5_new_dir, "lfc_plot.png"),
       device = "png", scale = 1.2, width = 2.5, height = 2, dpi = 330)

# plot discovery results
comparison_plot <- compare_calibration_and_discovery_results(calibration_result = calibration_result,
                                                             discovery_result = discovery_result)
comparison_plot <- comparison_plot + theme(legend.position = c(0.7, 0.1),
      legend.margin=margin(t = -0.5, unit='cm'),
      legend.title= element_blank(),) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5)))
  
volcano_plot <- make_volcano_plot(discovery_result = discovery_result,
                                  x_limits = c(-2.0, 2.0),
                                  transparency = 0.2, point_size = 0.5)
ggsave(filename = paste0(fig_5_new_dir, "comparison_plot.png"), plot = comparison_plot,
       device = "png", scale = 1.2, width = 2.5, height = 2, dpi = 330)
ggsave(filename = paste0(fig_5_new_dir, "volcano_plot.png"), plot = volcano_plot,
       device = "png", scale = 1.2, width = 2.5, height = 2, dpi = 330)