library(tidyverse)
library(katlabutils)
library(cowplot)
library(ondisc)

# load functions and data
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
pc_res <- readRDS(paste0(result_dir, "positive_control_analysis/pc_results_processed.rds"))
reject_thresh <- 1e-5
undercover_res <- readRDS(paste0(result_dir, "undercover_grna_analysis/undercover_result_grp_1_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF)
n_false_rejections <- undercover_res |>
  compute_n_bonf_rejections() |>
  rename("n_ntc_reject" = "n_reject") |>
  select(-Method)

##########
# TOP HALF
##########
n_reject_df <- pc_res |>
  filter(n_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_control >= N_NONZERO_CONTROL_CUTOFF) |>
  group_by(dataset, method) |>
  summarize(n_pc_reject = sum(p_value < reject_thresh),
            Method = Method[1]) |>
  group_by(dataset) |>
  mutate(n_pc_reject = ifelse(n_pc_reject == 0, 0.02 * max(n_pc_reject), n_pc_reject),
         Method = forcats::fct_recode(Method,
                                      "SCEPTRE" = "Sceptre",
                                      "Liscovitch Meth." = "Liscovitch Method",
                                      "Schraivogel Meth." = "Schraivogel Method",
                                      "Weissman Meth." = "Weissman Method"),
         Method = forcats::fct_relevel(Method,
                                       "SCEPTRE",
                                       "Seurat De",
                                       "Weissman Meth.",
                                       "Mimosca",
                                       "Liscovitch Meth.",
                                       "Schraivogel Meth.")) |>
  left_join(n_false_rejections,
            by = c("dataset", "method")) |>
  filter(n_ntc_reject < 50)

make_n_rejected_plot_for_dataset <- function(n_reject_df, dataset, tit, y_text = TRUE, extra_space = FALSE) {
  curr_n_reject_df <- n_reject_df |> filter(dataset == !!dataset)
  # max_rejected <- curr_n_reject_df$n_pc_reject |> max()
  
  integer_breaks <- function(n = 5, ...) {
    fxn <- function(x) {
      breaks <- floor(pretty(x, n, ...))
      names(breaks) <- attr(breaks, "labels")
      breaks
    }
    return(fxn)
  }
  
  p_0 <- curr_n_reject_df |>
    arrange(Method) |>
    ggplot2::ggplot(ggplot2::aes(x = Method,
                                 y = n_pc_reject,
                                 fill = Method)) +
    ggplot2::geom_col(col = "black") +
    ylab("N discoveries") +
    xlab("Method") +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values = my_cols) +
    scale_color_manual(values = my_cols) +
    scale_y_continuous(breaks = integer_breaks(),
                       expand = c(0, 0))
  
  legend <- get_legend(p_0)
    
  p <- p_0 +
    my_theme_no_legend +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(tit)
  
  if (!y_text) {
    p <- p + theme(axis.title.y = element_blank())
  }
    
  if (extra_space) {
    p <- p + theme(axis.title.y = element_text(margin = margin(t = 5, r = 10, b = 5, l = 0,  unit = "pt")))
  }
  
  return(list(p = p, legend = legend))
}

p_a <- make_n_rejected_plot_for_dataset(n_reject_df = n_reject_df, dataset = "frangieh_ifn_gamma_gene", tit = "Frangieh (IFN-\u03B3)", y_text = TRUE, extra_space = TRUE)
p_b <- make_n_rejected_plot_for_dataset(n_reject_df = n_reject_df, "frangieh_control_gene", "Frangieh (control)", y_text = TRUE, extra_space = TRUE)
p_c <- make_n_rejected_plot_for_dataset(n_reject_df = n_reject_df, "frangieh_co_culture_gene", "Frangieh (co culture)", extra_space = FALSE, y_text = FALSE)
p_d <- make_n_rejected_plot_for_dataset(n_reject_df = n_reject_df, "papalexi_eccite_screen_gene", "Papalexi (gene)",  extra_space = FALSE, y_text = FALSE)
p_e <- make_n_rejected_plot_for_dataset(n_reject_df = n_reject_df, "papalexi_eccite_screen_protein", "Papalexi (protein)",  extra_space = FALSE, y_text = FALSE)
p_f <- make_n_rejected_plot_for_dataset(n_reject_df = n_reject_df, "schraivogel_enhancer_screen", "Schraivogel",  extra_space = FALSE, y_text = FALSE)

fig_top <- plot_grid(p_a$p, p_c$p, p_e$p, p_b$p, p_d$p, p_f$p, nrow = 2, labels = c("a", "c", "e", "b", "d", "f"), align = "v")
legend <- p_e$legend
fig_top_2 <- plot_grid(fig_top, legend, ncol = 2, rel_widths = c(0.81, 0.19))

############
# PANELS g-h
############

make_p_val_vs_sample_size_plot <- function(pc_res_sub, tit, print_legend) {
    p <- pc_res_sub |>
    mutate(reject = ifelse(p_value < reject_thresh, "Signif.", "Not signif.")) |>
    ggplot(mapping = aes(x = n_treatment, y = p_value, col = reject)) +
    annotate("rect", xmin = -Inf, xmax = N_NONZERO_TREATMENT_CUTOFF, ymin = 0, ymax = Inf, fill = "slategray1") +
    annotate("rect", xmin = N_NONZERO_TREATMENT_CUTOFF, xmax = Inf, ymin = 0, ymax = Inf, fill = "lightpink") +
    geom_hline(yintercept = reject_thresh, linetype = "dashed", col = "darkred") +
    geom_point(alpha = 0.85, size = 0.95) +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 2),
                       breaks = c(0, 1, 10, 100), expand = c(0.02, 0)) +
    scale_y_continuous(trans = revlog_trans(10)) +
    xlab("N treatment cells with expression") +
    ylab("SCEPTRE p-value") +
    ggtitle(tit) +
    scale_color_manual(values = c("Signif." = "black", "Not signif." = "slategrey"),
                       breaks = c("Signif.", "Not signif."))
      
    if (print_legend) {
        p <- p +
          my_theme +
          theme(legend.title= element_blank(),
                legend.position = c(0.15, 0.87),
                legend.text=element_text(size = 9),
                legend.margin=margin(t = -0.5, b = 0, unit='cm'),
                legend.key = element_rect(fill = "slategray1"),
                legend.background = element_rect(fill = "slategray1")) +
          guides(color = guide_legend(
            keywidth = 0.0,
            keyheight = 0.15,
            default.unit="inch"))
          
      } else {
        p <- p + my_theme_no_legend
      }
}

p_g <- make_p_val_vs_sample_size_plot(pc_res |>
                                        filter(dataset == "frangieh_co_culture_gene",
                                               method == "sceptre"),
                                      "Frangieh (co culture) positive controls",
                                      print_legend = TRUE)

p_h <- make_p_val_vs_sample_size_plot(pc_res |>
                                         filter(dataset == "frangieh_control_gene", method == "sceptre"),
                                      "Frangieh (control) positive controls",
                                      print_legend = FALSE)
fig_bottom <- plot_grid(p_g, p_h, nrow = 1, labels = c("g", "h"))

############
# CREATE FIG
############
fig <- plot_grid(fig_top_2, fig_bottom, nrow = 2, rel_heights = c(0.6, 0.4), align = "v", axis = "l")

to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_5/fig_5.png")
ggsave(filename = to_save_fp, plot = fig, device = "png",
       scale = 1.1, width = 6.5, height = 6.25, dpi = 330)
