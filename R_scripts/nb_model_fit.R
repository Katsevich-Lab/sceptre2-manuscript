library(katlabutils)
library(tidyverse)
library(ondisc)
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
N_GENES <- 500
datasets <- c("frangieh/ifn_gamma",
              "frangieh/co_culture",
              "frangieh/control",
              "papalexi/eccite_screen")

res <- lapply(X = datasets, FUN = function(dataset) {
  mm <- lowmoi::load_dataset_multimodal(paper_fp = dataset, offsite_dir = sceptre2_dir)
  response_odm <- mm |> get_modality("gene")
  grna_assignments <- mm |>
    get_modality("grna_assignment") |>
    lowmoi:::get_target_assignments_via_max_op()
  response_odm_ntc <- response_odm[,grna_assignments == "non-targeting"]
  gene_ids <- response_odm_ntc |>
    get_feature_covariates() |>
    arrange(desc(n_nonzero)) |>
    slice(1:N_GENES) |>
    row.names()
  
  # perform the sample split
  set.seed(4)
  n_ntc_cells <- ncol(response_odm_ntc)
  samp_idxs <- sample(x = n_ntc_cells,
                      size = n_ntc_cells/2,
                      replace = FALSE)
  my_formula_str <- response_odm@misc$nb_regression_formula
  my_formula <- stats::as.formula(paste0("expression ", my_formula_str))
  cell_covariates_ntc <- response_odm_ntc |> get_cell_covariates()
  cell_covariates_ntc_s1 <- cell_covariates_ntc[samp_idxs,]
  cell_covariates_ntc_s2 <- cell_covariates_ntc[-samp_idxs,]
  
  thetas <- sapply(gene_ids, FUN = function(gene_id) {
    print(gene_id)
    expression <- as.numeric(response_odm_ntc[[gene_id, samp_idxs]])
    curr_data_matrix <- dplyr::mutate(cell_covariates_ntc_s1, expression = expression)
    theta <- lowmoi:::estimate_size(df = curr_data_matrix,
                                    formula = my_formula)
    theta <- min(max(theta, 0.1), 1000)
  })
  
  # p-values from deviance-based goodness of fit tests -- ideally the p-values are uniform
  fit_ps <- sapply(gene_ids, FUN = function(gene_id) {
    print(gene_id)
    curr_theta <- thetas[[gene_id]]
    expression <- as.numeric(response_odm_ntc[[gene_id, -samp_idxs]])
    curr_data_matrix <- dplyr::mutate(cell_covariates_ntc_s2, expression = expression)
    fit_nb <- glm(formula = my_formula,
                  family = MASS::neg.bin(curr_theta),
                  data = curr_data_matrix)
    fit_p <- pchisq(fit_nb$deviance,
                    df = fit_nb$df.residual,
                    lower.tail = FALSE)
  })
  data.frame(theta = thetas, p = fit_ps, dataset = dataset)
}) |> data.table::rbindlist()


res |>
  ggplot(mapping = aes(y = p, col = dataset)) +
  stat_qq_points(ymin = 1e-9, size = 0.8) +
  geom_abline(col = "darkred") +
  stat_qq_band() +
  theme_bw() +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "Expected quantile", y = "Observed quantile") +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
