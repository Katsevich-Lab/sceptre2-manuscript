pc_res <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                 "results/positive_control_analysis/pc_results_sceptre_unfiltered_0523_processed.rds") |>
  readRDS() |> na.omit()
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), 
                            "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
library(ggplot2)
library(dplyr)
library(katlabutils)
library(cowplot)
source(shared_fig_script)
reject_thresh <- 1e-5 
conflicts_prefer(dplyr::filter)
N_NONZERO_TREATMENT_CUTOFF <- 7

make_p_val_vs_sample_size_plot <- function(pc_res_sub, ylab = TRUE, xlab = TRUE) {
  tit <- as.character(pc_res_sub$dataset_rename[1])
  p <- pc_res_sub |>
    mutate(reject = ifelse(p_value < reject_thresh, "Signif.", "Not signif.")) |>
    mutate(p_value = ifelse(p_value < 1e-7, 1e-7, p_value)) |>
    ggplot(mapping = aes(x = n_treatment, y = p_value, col = reject)) +
    annotate("rect", xmin = -Inf, xmax = N_NONZERO_TREATMENT_CUTOFF, ymin = 0, ymax = Inf, fill = "slategray1") +
    annotate("rect", xmin = N_NONZERO_TREATMENT_CUTOFF, xmax = Inf, ymin = 0, ymax = Inf, fill = "lightpink") +
    geom_hline(yintercept = reject_thresh, linetype = "dashed", col = "darkred") +
    geom_point(alpha = 0.85, size = 0.95) +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 2),
                       breaks = c(0, 1, 10, 100), expand = c(0.02, 0)) +
    scale_y_continuous(trans = revlog_trans(10)) +
    xlab(if (xlab) "Effective sample size" else "") +
    ylab("SCEPTRE p-value") +
    ggtitle(tit) +
    scale_color_manual(values = c("Signif." = "black", "Not signif." = "slategrey"),
                       breaks = c("Signif.", "Not signif.")) +
    my_theme_no_legend +
    (if (!ylab) theme(axis.title.y = element_blank()))
  return(p)
}

p_a <- make_p_val_vs_sample_size_plot(pc_res |> filter(dataset == "frangieh_co_culture_gene"),
                                     ylab = TRUE, xlab = FALSE)
p_b <- make_p_val_vs_sample_size_plot(pc_res |> filter(dataset == "frangieh_control_gene"),
                                      ylab = FALSE, xlab = TRUE)
p_c <- make_p_val_vs_sample_size_plot(pc_res |> filter(dataset == "frangieh_ifn_gamma_gene"),
                                      ylab = FALSE, xlab = FALSE)

fig_bottom <- plot_grid(p_a, p_b, p_c, nrow = 1, labels = c("a", "b", "c"))

fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
             "sceptre2-manuscript/R_scripts/figure_creation/fig_s7/fig_s7.png")
ggsave(filename = fp, plot = fig_bottom, device = "png", scale = 1, width = 7, height = 3, dpi = 320)
