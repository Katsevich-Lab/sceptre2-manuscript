my_cols <- c("KS test" = "purple3",
             "MAST" = "aquamarine4", 
             "MIMOSCA" = "palevioletred3",
             "t-test" = "navy",
             "Seurat-Wilcox" = "dodgerblue3",
             "Seurat-NB" = "darkgoldenrod3",
             "Seurat-Wilcox (extreme filtering)" = "dodgerblue",
             "Seurat-Wilcox (standard filtering)" = "dodgerblue4",
             "NB regression (no covariates)" = "darkgoldenrod4",
             "NB regression (w/ covariates)" = "darkgoldenrod3",
             "SCEPTRE" = "firebrick3",
             "SCEPTRE (no covariates)" = "dodgerblue",
             "Permutation test" = "slategray4")

my_theme <- theme_bw() + theme(axis.line = element_line(color = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank(),
                               plot.title = element_text(hjust = 0.5, size=11))
my_theme_no_legend <- my_theme + theme(legend.position = "none")

compute_n_bonf_rejections <- function(undercover_res, alpha = 0.1) {
  n_bonf_rej <- undercover_res |>
    dplyr::group_by(dataset, method) |>
    dplyr::summarize(reject = (p_value < alpha/dplyr::n()),
                     Method = Method[1]) |>
    dplyr::summarize(n_reject = sum(reject),
                     Method = Method[1]) |>
    dplyr::ungroup()
  
  max_reject <- max(n_bonf_rej$n_reject)
  n_bonf_rej <- n_bonf_rej |>
    mutate(n_reject = ifelse(n_reject == 0, min(max_reject/50, 0.1), n_reject))
  
  return(n_bonf_rej)
}

N_NONZERO_TREATMENT_CUTOFF <- 7
N_NONZERO_CONTROL_CUTOFF <- 7

# define binning function
bin <- function(x, n_bins) {
  quantiles <- quantile(x = x, probs = seq(0, 1, 1/n_bins))
  stop <- quantiles[seq(2, length(quantiles))] - 1L
  start <- quantiles[seq(1, length(quantiles) - 1L)]
  labs <- paste0("[", start, ", ", stop, "]")
  out <- cut(x, breaks = quantile(x = x, probs = seq(0, 1, 1/n_bins)),
             include.lowest = TRUE, right = FALSE, labels = labs)
  return(out)
}
