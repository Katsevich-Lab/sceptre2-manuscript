
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
