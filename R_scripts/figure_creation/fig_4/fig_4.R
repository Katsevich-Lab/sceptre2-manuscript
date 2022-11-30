library(tidyverse)
library(katlabutils)
library(cowplot)
library(ondisc)

# load functions and data
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
undercover_res <- readRDS(paste0(result_dir,
                                 "undercover_grna_analysis/undercover_result_grp_1_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF,
         method %in% c("sceptre", "seurat_de")) |>
  mutate(Method = fct_recode(Method,
                             "SCEPTRE" = "Sceptre"))

########################################################
# PANNELS A-D: SCEPTRE vs. Seurat De on several datasets
########################################################
my_values <- my_cols[names(my_cols) %in% c("Seurat De", "SCEPTRE")]
get_plots_for_dataset <- function(df_sub, tit, print_legend, legend_position = c(0.45, 0.85)) {
p_qq <- ggplot(data = df_sub, mapping = aes(y = p_value, col = Method)) +
    stat_qq_points(ymin = 1e-8, size = 0.85) +
    stat_qq_band() +
    scale_x_continuous(trans = revlog_trans(10)) +
    scale_y_continuous(trans = revlog_trans(10)) +
    labs(x = "Expected null p-value", y = "Observed p-value") +
    geom_abline(col = "black") +
    ggtitle(tit) +
    scale_color_manual(values = my_values,
                       drop = FALSE,
                       breaks = c("SCEPTRE",
                                  "Seurat De",
                                  "SCEPTRE (w/o covariates)",
                                  "NB Regression"))
  
  if (print_legend) {
    p_qq <- p_qq +
      my_theme +
      theme(legend.title= element_blank(),
            legend.position = legend_position,
            legend.text=element_text(size = 11),
            legend.margin=margin(t = 0, unit='cm')) +
      guides(color = guide_legend(
        keywidth = 0.0,
        keyheight = 0.2,
        default.unit="inch"))
  } else {
    p_qq <- p_qq + my_theme_no_legend
  }
  
  n_bonf_rej <- df_sub |> compute_n_bonf_rejections()
  max_reject <- max(n_bonf_rej$n_reject)
  
  breaks_v <-  seq(0, max_reject, by = if (max_reject >= 7) 2 else 1)
  p_bar <- n_bonf_rej |> ggplot2::ggplot(ggplot2::aes(x = Method, y = n_reject, fill = Method)) +
    ggplot2::geom_col(col = "black") +
    ylab("N Bonferoni rejections") +
    xlab("Method") + my_theme_no_legend +
    theme(axis.text.x = element_blank(),
          plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt")) +
    scale_y_continuous(breaks = breaks_v, expand = c(0, 0)) +
    ggtitle("") +
    scale_fill_manual(values = my_values)
  
  return(list(p_qq = p_qq, p_bar = p_bar))
}

# 1.
papa_plots <- get_plots_for_dataset(undercover_res |> filter(dataset == "papalexi_eccite_screen_gene"),
                                    "Papalexi (gene) neg. controls",
                                    print_legend = FALSE)
# 2.
papa_protein_plots <- get_plots_for_dataset(undercover_res |>
                                              filter(dataset == "papalexi_eccite_screen_protein",
                                                     method %in% c("sceptre", "seurat_de")),
                                            "Papalexi (protein) neg. controls",
                                            print_legend = FALSE)
# 3.
ifn_gama_plots <- get_plots_for_dataset(undercover_res |>
                                          filter(dataset == "frangieh_ifn_gamma_gene"),
                                        "Frangieh (IFN-\u03B3) neg. controls",
                                        print_legend = TRUE,
                                        legend_position = c(0.25, 0.8))
# 4.
co_culture_plots <- get_plots_for_dataset(undercover_res |>
                                            filter(dataset == "frangieh_co_culture_gene",
                                                   method %in% c("sceptre", "seurat_de")),
                                            "Frangieh (co culture) neg. controls",
                                          print_legend = FALSE)

# 5.
control_plots <- get_plots_for_dataset(undercover_res |> 
                                         filter(dataset == "frangieh_control_gene",
                                                method %in% c("sceptre", "seurat_de")),
                                       "Frangieh (control) neg. controls",
                                       print_legend = FALSE)
# 6
enh_screen <- get_plots_for_dataset(undercover_res |> 
                                      filter(dataset == "schraivogel_enhancer_screen",
                                             method %in% c("sceptre", "seurat_de")),
                                    "Schraivogel neg. controls",
                                    print_legend = FALSE)

fig <- cowplot::plot_grid(ifn_gama_plots$p_qq, ifn_gama_plots$p_bar,
                          co_culture_plots$p_qq, co_culture_plots$p_bar,
                          control_plots$p_qq, control_plots$p_bar,
                          papa_plots$p_qq, papa_plots$p_bar,
                          papa_protein_plots$p_qq, papa_protein_plots$p_bar,
                          enh_screen$p_qq, enh_screen$p_bar,
                          labels = c("a", "", "b", "", "c", "", "d", "", "e", "", "f", ""),
                          rel_widths = c(0.38, 0.12, 0.38, 0.12),
                          ncol = 4, nrow = 3, align = "h")

to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_3/r_output.png")
ggsave(filename = to_save_fp, plot = fig, device = "png",
       scale = 1.1, width = 7, height = 7, dpi = 330)
