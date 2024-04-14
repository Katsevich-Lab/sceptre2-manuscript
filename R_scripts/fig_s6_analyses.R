library(katlabutils)
library(tidyverse)
library(ondisc)
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
N_GENES <- 1000

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
