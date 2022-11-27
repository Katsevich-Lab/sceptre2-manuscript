library(katlabutils)
library(tidyverse)
library(ondisc)
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
N_GENES <- 1000

#########################################
# GOODNESS OF FIT TESTS FOR NB REGRESSION
#########################################
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
    if (mean(expression >= 3) >= 0.95) {
      curr_data_matrix <- dplyr::mutate(cell_covariates_ntc_s1, expression = expression)
      theta <- lowmoi:::estimate_size(df = curr_data_matrix,
                                      formula = my_formula)
      theta <- min(max(theta, 0.1), 1000) 
    } else {
      NA
    }
  })
  
  fit_ps <- sapply(gene_ids, FUN = function(gene_id) {
    print(gene_id)
    curr_theta <- thetas[[gene_id]]
    if (!is.na(curr_theta)) {
      expression <- as.numeric(response_odm_ntc[[gene_id, -samp_idxs]])
      curr_data_matrix <- dplyr::mutate(cell_covariates_ntc_s2, expression = expression)
      fit_nb <- glm(formula = my_formula,
                    family = MASS::neg.bin(curr_theta),
                    data = curr_data_matrix)
      fit_p <- pchisq(fit_nb$deviance,
                      df = fit_nb$df.residual,
                      lower.tail = FALSE)
    } else {
      fit_p <- NA
    }
    return(fit_p)
  })
  data.frame(theta = thetas, p = fit_ps,
             dataset = dataset, response_id = gene_ids)
}) |> data.table::rbindlist()

res <- res |>
  lowmoi:::replace_slash_w_underscore() |>
  na.omit() |>
  mutate(dataset = paste0(dataset, "_gene")) |>
  as.data.frame()
saveRDS(object = res,
        file = paste0(result_dir, "extra_analyses/goodness_of_fit_tests.rds"))

################################################################################
# Testing for association between gRNA presence/absence and biological replicate
################################################################################
# load papalexi data
response_odm <- lowmoi::load_dataset_modality("papalexi/eccite_screen/gene")
grna_odm <- lowmoi::load_dataset_modality("papalexi/eccite_screen/grna_assignment")

# find a gRNA group that is strongly associated with bio_rep
grna_assignments <- lowmoi:::get_grna_assignments_via_max_op(grna_odm)
nt_cells <- grepl(pattern = "^NTg*", x = grna_assignments)
nt_grna_assignments <- grna_assignments[nt_cells]
unique_nt_grnas <- unique(nt_grna_assignments)
biorep_vect <- response_odm |>
  get_cell_covariates() |>
  pull(bio_rep)
biorep_vect_nt <- biorep_vect[nt_cells]
fisher_exact_p <- sapply(unique_nt_grnas, function(unique_nt_grna) {
  grna_binary_vect <- as.integer(unique_nt_grna == nt_grna_assignments)
  cont_table <- as.matrix(table(biorep_vect_nt, grna_binary_vect))
  fisher.test(x = cont_table)$p.value
})
saveRDS(object = fisher_exact_p,
        file = paste0(result_dir, "extra_analyses/papalexi_grna_confounding_tests.rds"))

#####################################################################################
# Testing for association between (relative) gene expression and biological replicate
#####################################################################################
gene_exp_mat <- as.matrix(response_odm[[,nt_cells]])
rownames(gene_exp_mat) <- get_feature_ids(response_odm)
cell_cov <- (response_odm |> get_cell_covariates())[nt_cells,]
gene_ids <- response_odm |>
  get_feature_covariates() |>
  arrange(desc(mean_expression)) |>
  slice(seq(1, N_GENES)) |>
  row.names()
full_formula <- formula(expressions ~ bio_rep + log(n_umis))
reduced_formula <-  formula(expressions ~ log(n_umis))

lrt_p <- sapply(X = gene_ids, function(gene_id) {
  print(paste0("Fitting model for ", gene_id))
  expressions <- gene_exp_mat[gene_id,]
  curr_cell_cov <- mutate(cell_cov, expressions = expressions)
  # estimate the size
  est_size <- lowmoi:::estimate_size(df = curr_cell_cov,
                                     formula = full_formula)
  # fit the NB regression models
  full_nb_reg <- glm(formula = full_formula,
                     family = MASS::negative.binomial(est_size),
                     data = curr_cell_cov)
  reduced_nb_reg <- glm(formula = reduced_formula,
                        family = MASS::negative.binomial(est_size),
                        data = curr_cell_cov)
  fit <- anova(reduced_nb_reg, full_nb_reg, test = "LRT")
  fit$`Pr(>Chi)`[2]
}) # All LRTs are highly significant, indicating a strong association between (relative) gene expression and biological replicate. I will (somewhat arbitrarily) choose the second most highly expressed gene to plot.
saveRDS(object = lrt_p, file = paste0(result_dir, "extra_analyses/papalex_gene_confounding_tests.rds"))
