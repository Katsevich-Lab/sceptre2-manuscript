my_cols <- c("Weissman Method" = "purple3",
             "Weissman Meth." = "purple3",
             "Schraivogel Method" = "dodgerblue3",
             "Schraivogel Meth." = "dodgerblue3",
             "Mimosca" = "orange3",
             "Liscovitch Method" = "darkslategray4",
             "Liscovitch Meth." = "darkslategray4",
             "Seurat De" = "lightseagreen",
             "Seurat De (w/ strict QC)" = "dodgerblue4",
             "Seurat De (w/o strict QC)" = "dodgerblue",
             "NB Reg (w/ strict QC)" = "indianred4",
             "NB Reg (w/o strict QC)" = "indianred1",
             "NB Reg (w/o covariates)" = "indianred4",
             "NB Reg (w/ covariates)" = "indianred1",
             "NB Regression" = "dodgerblue3",
             "SCEPTRE" = "firebrick3",
             "SCEPTRE (w/o covariates)" = "orange3",
             "Permutation test" = "darkslategray4")

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