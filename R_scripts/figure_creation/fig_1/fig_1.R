# Load packages
library(tidyverse)
library(katlabutils)
library(cowplot)

# Load scripts and results
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
undercover_res <- readRDS(paste0(result_dir,
                                 "undercover_grna_analysis/undercover_result_grp_1_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF)

# The following plot (the letters are rows):
# a) three empty columns for undercover schematic
# b) col1: Frangieh IFN gamma transformed
#    col2: Frangieh IFN gamma transformed
#    col3: Frangieh IFN gamma N bonferoni rejections (at level 0.1)
# c) same as above, but with Papalexi
# d) same as above, but with simulated OR Schraivogel enhancer screen
# restricting attention to pairs with >= 10 treatment cells and > 10 control cells in all cases

make_figure_row <- function(dataset, name, print_legend) {
  my_methods <- c("Weissman Method", "Schraivogel Method", "Mimosca", "Liscovitch Method", "Seurat De")
  my_values <- my_cols[names(my_cols) %in% my_methods]
  
  df_sub <- undercover_res |>
    filter(dataset == !!dataset,
           Method %in% my_methods)
  
  p1 <- ggplot(data = df_sub, mapping = aes(y = p_value, col = Method)) +
    stat_qq_points(ymin = 1e-8, size = 0.55) +
    stat_qq_band() +
    scale_x_reverse() +
    scale_y_reverse() +
    labs(x = "Expected null p-value", y = "Observed p-value") +
    geom_abline(col = "black") +
    scale_color_manual(values = my_values) +
    ggtitle(paste0("QQ-plot (", name, ")"))
  
  if (print_legend) {
    p1 <- p1 +
      my_theme +
      theme(legend.title = element_blank(),
            legend.position = c(0.71, 0.2),
            legend.text=element_text(size = 9),
            legend.margin=margin(t = 0, unit='cm')) +
      guides(color = guide_legend(
        keywidth = 0.0,
        keyheight = 0.1,
        default.unit="inch",
        override.aes = list(size = 2.5)))
  } else {
    p1 <- p1 + my_theme_no_legend
  }
  
  p2 <- ggplot(data = df_sub, mapping = aes(y = p_value, col = Method)) +
    stat_qq_points(ymin = 1e-8, size = 0.55) +
    stat_qq_band() +
    scale_x_continuous(trans = revlog_trans(10)) +
    scale_y_continuous(trans = revlog_trans(10)) +
    labs(x = "Expected null p-value", y = "Observed p-value") +
    geom_abline(col = "black") + my_theme_no_legend +
    scale_color_manual(values = my_values) +
    ggtitle("Transformed QQ-plot")
  
  bonf_reject_df <- compute_n_bonf_rejections(df_sub)
  
  p3 <- bonf_reject_df |>
    ggplot2::ggplot(ggplot2::aes(x = Method, y = n_reject, fill = Method)) +
    ggplot2::geom_col(col = "black") +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.2),
                       expand = c(0,0),
                       breaks = c(0, 1, 10, 100, 1000, 8000)) +
    # scales::pseudo_log_trans() +
    # scale_y_log10(expand = c(0, 0)) +
    ylab("N Bonferoni rejections") +
    xlab("Method") + my_theme_no_legend +
    theme(axis.text.x = element_blank()) +
    scale_fill_manual(values = my_cols) +
    ggtitle("N rejections")
  
  p_row <- plot_grid(p1, p2, p3, nrow = 1, labels = NULL,
                     align = "h", rel_widths = c(0.42, 0.42, 0.26))
  p_row
}

r0 <- ggplot() +
  theme_minimal() +
  ggtitle("Undercover gRNA calibration assessment") +
  theme(plot.title = element_text(hjust = 0.5, size=11))
r1 <- make_figure_row(dataset = "papalexi_eccite_screen_gene", name = "Papalexi gene modality", print_legend = TRUE)
r2 <- make_figure_row(dataset = "frangieh_ifn_gamma_gene", name = "Frangieh IFN-\u03B3", print_legend = FALSE)

fig <- plot_grid(r0, r1, r2, nrow = 3,
                 labels = c("a", "b", "c"),
                 rel_heights = c(0.3, 0.35, 0.35))

to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_1/r_output.png")
ggsave(filename = to_save_fp, plot = fig, device = "png",
       scale = 1.1, width = 6.5, height = 6.75, dpi = 330)
