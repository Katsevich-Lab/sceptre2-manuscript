my_cols <- c("Weissman Method" = "purple3",
             "Schraivogel Method" = "dodgerblue3",
             "Mimosca" = "orange3",
             "Liscovitch Method" = "darkslategray4",
             "Seurat De" = "lightseagreen",
             "Seurat De (w/ strict QC)" = "dodgerblue4",
             "Seurat De (w/o strict QC)" = "dodgerblue",
             "NB Reg (w/ strict QC)" = "indianred4",
             "NB Reg (w/o strict QC)" = "indianred1",
             "NB Reg (w/o covariates)" = "indianred4",
             "NB Reg (w/ covariates)" = "indianred1",
             "NB Regression" = "dodgerblue3",
             "SCEPTRE" = "firebrick3",
             "SCEPTRE (w/o covariates)" = "orange3")

my_theme <- theme_bw() + theme(axis.line = element_line(color = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank(),
                               plot.title = element_text(hjust = 0.5, size=11))
my_theme_no_legend <- my_theme + theme(legend.position = "none")

compute_n_bonf_rejections <- function(undercover_res, alpha = 0.1) {
  out <- undercover_res |>
    dplyr::group_by(dataset_rename, Method) |>
    dplyr::summarize(reject = (p_value < alpha/dplyr::n())) |>
    dplyr::summarize(n_reject = sum(reject)) |>
    dplyr::ungroup()
  return(out)
}
