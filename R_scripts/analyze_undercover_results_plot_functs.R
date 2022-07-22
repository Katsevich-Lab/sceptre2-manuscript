# define several functions
# 1. combine schraivogel enhancer screens
combine_schraivogel_enhancer_screens <- function(undercover_res) {
  undercover_res |> dplyr::mutate(dataset = dataset |> forcats::fct_recode(schraivogel_enhancer_screen = "schraivogel_enhancer_screen_chr11_gene",
                                                                           schraivogel_enhancer_screen = "schraivogel_enhancer_screen_chr8_gene"))
}

# 2. perform sanity check
perform_sanity_check <- function(undercover_res) {
  # i) confirm number of p-values coincides across methods for a given dataset
  n_p_vals_coincide <- undercover_res |>
    dplyr::group_by(dataset, method) |>
    dplyr::summarize(count = dplyr::n()) |>
    dplyr::summarize(n_pvals_coincide = all(diff(count) == 0)) |>
    dplyr::pull(n_pvals_coincide) |> all()
  
  # i) confirm there are no NAs
  no_nas <- all(undercover_res |>
        dplyr::group_by(dataset, method) |>
        dplyr::summarize(count = sum(is.na(p_value))) |>
        dplyr::pull(count) == 0)
  return(n_p_vals_coincide & no_nas)  
}

# 3. update dataset names
update_dataset_names <- function(undercover_res) {
  out <- undercover_res |>
    dplyr::group_by(dataset, method) |>
    dplyr::mutate(n_pairs = dplyr::n(),
                  dataset_rename = stringr::str_to_title(gsub(pattern = "_",replacement = " ", x = dataset)),
                  dataset_rename_w_pairs = paste0(dataset_rename, " (", n_pairs[1], " pairs)"),
                  Method = stringr::str_to_title(gsub(pattern = "_",replacement = " ", x = method)))
  return(out)
}
  
# 4. create the untransformed qq-plot
make_trans_qq_plot <- function(undercover_res) {
  p <- ggplot(data = undercover_res, mapping = aes(y = p_value, color = Method)) +
    stat_qq_points(ymin = 1e-10, size = 0.8) +
    facet_wrap(~dataset_rename_w_pairs, scales = "free", labeller = label_wrap_gen(35)) +
    geom_abline() +
    stat_qq_band() +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    theme_bw() +
    labs(x = "Expected quantile", y = "Observed quantile") +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  return(p)
}

# 5. make transformed qq-plot
make_untrans_qq_plot <- function(undercover_res) {
  p <- ggplot(data = undercover_res, mapping = aes(y = p_value, color = Method)) +
    geom_vline(xintercept = 0.01) +
    stat_qq_points(ymin = 1e-10, size = 0.8) +
    facet_wrap(~dataset_rename_w_pairs, scales = "free", labeller = label_wrap_gen(35)) +
    geom_abline() +
    stat_qq_band() +
    scale_x_reverse() +
    scale_y_reverse() +
    theme_bw() +
    labs(x = "Expected quantile", y = "Observed quantile") +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  return(p)
}


make_computation_plot <- function(undercover_res) {
  comp_df <- undercover_res |>
    dplyr::select(undercover_grna, Dataset = dataset_rename, Method, clock_time, max_ram) |>
    dplyr::distinct() |>
    dplyr::group_by(Dataset, Method) |>
    dplyr::summarize(m_clock_time = mean(clock_time)/60, m_max_ram = mean(max_ram)) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(cols = c("m_clock_time", "m_max_ram"),
                        names_to = "metric", values_to = "value") |>
    dplyr::mutate(metric = factor(x = metric, levels = c("m_clock_time", "m_max_ram"),
                                  labels = c("Time (m)", "RAM (GB)"))) |>
    dplyr::mutate(value = ifelse(value < 1, 1.05, value))
  my_fct_order <- comp_df |>
    dplyr::filter(Method == "Schraivogel Method", metric == "RAM (GB)") |>
    dplyr::arrange(value) |>
    dplyr::pull(Dataset)
  comp_df <- comp_df |>
    dplyr::mutate(Dataset = factor(x = Dataset, levels = my_fct_order, labels = my_fct_order)) |>
    dplyr::filter(Dataset != "Liscovitch Method")
  
  p_undercover_comp <- ggplot(data = comp_df, mapping = aes(x = Method, y = value, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge", col = "black") +
    facet_wrap(.~metric, scales = "free_y") +
    ylab("") + xlab("Method") +
    theme_bw() +
    scale_y_continuous(trans='log2')
  return(p_undercover_comp)
}


# other error metrics: number of rejected pairs after Bonf correction
if (FALSE) {
  # save plots
  ggsave(filename = paste0(fig_dir, "undercov_grna_trans_grp3.png"),
         plot = p_undercover_stat_trans, device = "png", scale = 1, width = 11, height = 6, dpi = 330)
  ggsave(filename = paste0(fig_dir, "undercov_grna_untrans_grp3.png"),
         plot = p_undercover_stat_untrans, device = "png", scale = 1, width = 11, height = 6, dpi = 330)
  ggsave(filename = paste0(fig_dir, "undercov_grna_comp.pdf"),
         plot = p_undercover_comp, device = "pdf", scale = 0.8, width = 14, height = 5, dpi = 330)
  
  
  alpha <- 0.1
  bonf_correct <- res_proc |>
    dplyr::group_by(dataset, method) |>
    dplyr::summarize(reject = (p_value < alpha/dplyr::n())) |>
    dplyr::summarize(n_reject = sum(reject)) |>
    dplyr::ungroup()
  
  bonf_correct |>
    ggplot2::ggplot(ggplot2::aes(x = method, y = n_reject)) +
    ggplot2::geom_col() + 
    facet_wrap(~dataset) +
    ggplot2::scale_y_log10() +
    theme_bw()
}