my_cols <- c("Weissman Method" = "purple3",
             "Weissman Meth." = "purple3",
             "Schraivogel Method" = "aquamarine4", 
             "Schraivogel Meth." = "aquamarine4",
             "Mimosca" = "palevioletred3",
             "Liscovitch Method" = "navy",
             "Liscovitch Meth." = "navy",
             "Seurat De" = "dodgerblue3",
             "Seurat De (w/ strict QC)" = "dodgerblue4",
             "Seurat De (w/o strict QC)" = "dodgerblue",
             "NB Reg (w/o covariates)" = "darkgoldenrod3",
             "NB Reg (w/ covariates)" = "darkgoldenrod4",
             "NB Regression" = "goldenrod3",
             "SCEPTRE" = "firebrick3",
             "SCEPTRE (w/o covariates)" = "dodgerblue",
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