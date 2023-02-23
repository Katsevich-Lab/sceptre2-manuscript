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
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF) |>
  mutate(Method = forcats::fct_recode(Method,
                                      "SCEPTRE" = "Sceptre",
                                      "t-test" = "Liscovitch Method",
                                      "MAST" = "Schraivogel Method",
                                      "KS test" = "Weissman Method",
                                      "MIMOSCA" = "Mimosca"
  ))

alpha <- 0.1

my_cols <- c("Weissman Method" = "purple3",
             "Weissman Meth." = "purple3",
             "KS test" = "purple3",
             "Schraivogel Method" = "aquamarine4", 
             "Schraivogel Meth." = "aquamarine4",
             "MAST" = "aquamarine4",
             "Mimosca" = "palevioletred3",
             "MIMOSCA" = "palevioletred3",
             "Liscovitch Method" = "navy",
             "Liscovitch Meth." = "navy",
             "t-test" = "navy",
             "Seurat De" = "dodgerblue3",
             "Seurat De (w/ strict QC)" = "dodgerblue4",
             "Seurat De (w/o strict QC)" = "dodgerblue",
             "NB Reg (w/o covariates)" = "darkgoldenrod3",
             "NB Reg (w/ covariates)" = "darkgoldenrod4",
             "NB Regression" = "goldenrod3",
             "SCEPTRE" = "firebrick3",
             "SCEPTRE (w/o covariates)" = "dodgerblue",
             "Permutation test" = "slategray4")

# The following plot (the letters are rows):
# a) three empty columns for undercover schematic
# b) col1: Frangieh IFN gamma transformed
#    col2: Frangieh IFN gamma transformed
#    col3: Frangieh IFN gamma N bonferoni rejections (at level 0.1)
# c) same as above, but with Papalexi
# d) same as above, but with simulated OR Schraivogel enhancer screen
# restricting attention to pairs with >= 10 treatment cells and > 10 control cells in all cases

make_figure_row <- function(dataset, name, print_legend) {
  my_methods <- c("KS test", "MAST", "MIMOSCA", "t-test", "Seurat De")
  my_values <- my_cols[names(my_cols) %in% my_methods]
  
  df_sub <- undercover_res |>
    filter(dataset == !!dataset,
           Method %in% my_methods)
  
  bonf_thresh <- alpha/(nrow(df_sub)/length(my_methods))
  
  p2 <- ggplot(data = df_sub, mapping = aes(y = p_value, col = Method)) +
    stat_qq_points(ymin = 1e-8, size = 0.55) +
    stat_qq_band() +
    scale_x_continuous(trans = revlog_trans(10)) +
    scale_y_continuous(trans = revlog_trans(10)) +
    labs(x = "Expected null p-value", y = "Observed p-value") +
    geom_abline(col = "black") + 
    geom_hline(yintercept = bonf_thresh, linetype = "dashed") +
    scale_color_manual(values = my_values)
  
  if (print_legend) {
    p2 <- p2 +
      my_theme +
      theme(legend.title = element_blank(),
            legend.position = c(0.78, 0.22),
            legend.text=element_text(size = 9),
            legend.margin=margin(t = 0, unit='cm'),
            legend.background = element_blank()) +
      guides(color = guide_legend(
        keywidth = 0.0,
        keyheight = 0.1,
        default.unit="inch",
        override.aes = list(size = 2.5)))
  } else {
    p2 <- p2 + my_theme_no_legend
  }
  
  bonf_reject_df <- compute_n_bonf_rejections(df_sub, alpha = alpha)
  
  p3 <- bonf_reject_df |>
    ggplot2::ggplot(ggplot2::aes(x = Method, y = n_reject, fill = Method)) +
    ggplot2::geom_col(col = "black") +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.2),
                       expand = c(0,0),
                       breaks = c(0, 1, 10, 100, 1000, 8000)) +
    ylab("False positives") +
    xlab("Method") + my_theme_no_legend +
    theme(axis.text.x = element_blank()) +
    scale_fill_manual(values = my_cols)
  
  p_row <- plot_grid(p2, p3, nrow = 1, labels = NULL,
                     align = "h", rel_widths = c(2,1))
  p_row
}

r0 <- ggplot() +
  theme_minimal() +
  ggtitle("Undercover gRNA calibration assessment") +
  theme(plot.title = element_text(hjust = 0.5, size=11))

r1 <- make_figure_row(dataset = "papalexi_eccite_screen_gene", name = "Papalexi gene modality", print_legend = TRUE)
r2 <- make_figure_row(dataset = "frangieh_ifn_gamma_gene", name = "Frangieh IFN-\u03B3", print_legend = FALSE)

fig <- plot_grid(r1, r2, nrow = 1)


to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/genes_figs/fig_1_bottom.png")
ggsave(filename = to_save_fp, plot = fig, device = "png",
       scale = 1.1, width = 7, height = 2, dpi = 330)