my_cols <- ggsci::pal_npg()(7)

# define several functions

# 0. replace underscore with slash in dataset column
replace_dataset_underscore_with_slash <- function(datasets) {
  sub(pattern = "_", replacement = "/", x = datasets, fixed = TRUE) |>
    stringi::stri_replace_last_fixed(str = _, pattern = "_", replacement = "/")
}

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
update_dataset_names <- function(undercover_res, add_n_pairs = TRUE) {
  out <- undercover_res |>
    dplyr::group_by(dataset, method) |>
    dplyr::mutate(dataset_rename = stringr::str_to_title(gsub(pattern = "_",replacement = " ", x = dataset)),
                  Method = stringr::str_to_title(gsub(pattern = "_", replacement = " ", x = method)))
  if (add_n_pairs) {
    out <- out |> dplyr::mutate(n_pairs = dplyr::n(),
                                dataset_rename_w_pairs = paste0(dataset_rename, " (", n_pairs[1], " pairs)"))
  }
  out <- out |> dplyr::ungroup() |> dplyr::mutate(dataset_rename = NULL)
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

# 6. make computation plot
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

# 7. function to compute n Bonferoni-corrected hypotheses that have been rejected
compute_n_bonf_rejected <- function(undercover_res, alpha = 0.05) {
  # set.seed(10)
  # undercover_res_sub <- undercover_res |>
  #  dplyr::group_by(dataset_rename, Method) |>
  #  dplyr::mutate(n_hyp = dplyr::n()) |>
  #  dplyr::filter(n_hyp >= 5000) |>
  #  dplyr::sample_n(size = 5000)
  out <- undercover_res |>
    dplyr::group_by(dataset_rename_w_pairs, Method) |>
    dplyr::summarize(reject = (p_value < alpha/dplyr::n())) |>
    dplyr::summarize(n_reject = sum(reject)) |>
    dplyr::ungroup()
  return(out)
}


# 8. make n rejected pairs plot
make_n_rejected_pairs_plot <- function(n_rejected_df, y_max = 1e5, scales = "fixed", log_trans = TRUE) {
  n_rejected_df |>
    ggplot2::ggplot(ggplot2::aes(x = Method, y = n_reject, fill = Method)) +
    ggplot2::geom_col(col = "black") + 
    facet_wrap(~dataset_rename_w_pairs, labeller = label_wrap_gen(35), scales = scales) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title= element_blank()) +
    ylab("N rejected (after Bonf. correction)") +
    if (log_trans) ggplot2::scale_y_log10(limits = c(1, y_max)) else NULL 
}


# 9. make histograms
make_histograms <- function(undercover_res_to_plot, fill) {
  ggplot(data = undercover_res_to_plot, mapping = aes(x = p_value)) +
    facet_wrap(dataset_rename_w_pairs ~ ., scales = "free_y", labeller = label_wrap_gen(35)) +
    geom_histogram(bins = 15, aes(y = stat(density)),
                   col = "black", boundary = 0, fill = fill) +
    theme_bw() + xlab("p-value") + ylab("Density")
}


# 10. append the number of undercover cells
append_n_undercover_cells <- function(undercover_res) {
  undercover_res |>
    group_by(undercover_grna, dataset_slash) |>
    group_modify(function(tibble, key) {
      print(paste0("Working on ", key$undercover_grna, " ", key$dataset_slash))
      dataset_name <- as.character(key$dataset_slash)
      undercover_grna <- strsplit(x = as.character(key$undercover_grna), split = ",", fixed = TRUE) |> unlist()
      grna_feature_covariates <- lowmoi::get_grna_dataset_name(dataset_name, "assignment") |> 
        lowmoi::load_dataset_modality() |>
        ondisc::get_feature_covariates()
      n_undercover_cells <- grna_feature_covariates[undercover_grna, "n_nonzero"] |> sum()
      mutate(tibble, n_undercover_cells = n_undercover_cells)
    })
}
  
# 11. for each dataset-method pair, return the number of rejections and number of cells for each undercover gRNA
get_n_reject_n_undercover_cells_df <- function(undercover_res, alpha) {
  undercover_res |>
    dplyr::group_by(dataset_rename_w_pairs, method) |>
    dplyr::mutate(reject = (p_value < alpha/dplyr::n())) |>
    dplyr::group_by(dataset_rename_w_pairs, method, undercover_grna) |>
    dplyr::summarize(n_reject = sum(reject),
                     n_undercover_cells = n_undercover_cells[1],
                     dataset = dataset[1]) |>
    dplyr::ungroup()
}

# 12. associate n undercover cells w n rejections
associate_n_undercover_cells_w_n_rejections <- function(n_rejected_n_cells_df_sub) {
  lm_df <- n_rejected_n_cells_df_sub |>
    group_by(dataset_rename_w_pairs, method) |>
    group_modify(function(tibble, key) {
      fit <- lm(formula = n_reject ~ n_undercover_cells, data = tibble)
      s <- summary(fit)
      p <- s$coefficients["n_undercover_cells", "Pr(>|t|)"]
      coefs <- coef(fit)
      data.frame(intercept = coefs[[1]], slope = coefs[[2]], p_val = p)
    })
  
  p <- ggplot(data = n_rejected_n_cells_df_sub,
         mapping = aes(y = n_reject, x = n_undercover_cells)) +
    facet_wrap(. ~ dataset_rename_w_pairs, scales = "free", labeller = label_wrap_gen(35)) +
    geom_point() +
    geom_abline(data = lm_df,
                mapping = aes(slope = slope, intercept = intercept, col = p_val < 0.05),
                lwd = 0.7) + theme_bw() + 
    xlab("N undercover cells") + ylab("N reject")
  return(p)
}

# 13. plot gRNA-specific qq-plots, coloring the points according to whether they deviate from the expected uniform distribution according to a KS test
make_grna_specific_plots <- function(undercover_res_sub) {
  # perform gRNA-wise KS test; uniform tail cdf function
  tail_cdf <- function(x) {
    sapply(X = x, FUN = function(curr_x) {
      if (curr_x < 0) 0 else if (curr_x > 0.1) 1 else 10 * curr_x
    })
  }
  undercover_w_ks <- undercover_res_sub |>
    group_by(undercover_grna) |>
    group_modify(.f = function(tibble, key) {
      p_vals <- tibble$p_value
      p_vals <- p_vals[p_vals <= 0.1]
      ks_fit <- suppressWarnings(ks.test(p_vals, tail_cdf))
      ks_p <- ks_fit$p.value
      ks_d <- ks_fit$statistic
      mutate(tibble, ks_p_val = ks_p, ks_d_val = ks_d)
    }) |>
    mutate(ks_test = factor(ks_p_val < 0.001, levels = c(TRUE, FALSE),
                            labels = c("non-unif", "unif")))
  # create the plot
  p <- ggplot(data = undercover_w_ks, mapping = aes(y = p_value)) +
    stat_qq_points(ymin = 1e-10, size = 0.8, mapping = aes(col = ks_test)) +
    facet_wrap(~undercover_grna, scales = "free") +
    geom_abline() +
    stat_qq_band() +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    theme_bw() +
    labs(x = "Expected quantile", y = "Observed quantile") +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_color_manual(values = c("darkred", "darkblue"))
  return(p)
}
