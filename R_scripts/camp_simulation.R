library(ondisc)
library(sceptre2)
library(microbenchmark)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/extra_analyses/")

###########################
# STEP 1: SET UP SIMULATION
###########################

# 1. load Papalexi data
response_odm <- lowmoi::load_dataset_modality("papalexi/eccite_screen/gene")
grna_odm <- lowmoi::load_dataset_modality("papalexi/eccite_screen/grna_assignment")
grna_targets <- lowmoi::get_target_assignments_via_max_op(grna_odm)

# 2. get the gRNA group info
my_grna <- "CUL3"
grna_group_info <- sceptre2:::get_grna_group_info(grna_group_assignments = grna_targets,
                                                  input_grna_groups = my_grna)
idxs <- c(grna_group_info$grna_specific_idxs[[my_grna]],
          grna_group_info$grna_specific_idxs[["non-targeting"]])
orig_x <- c(rep(1, grna_group_info$n_cells_per_grna[[my_grna]]),
            rep(0, grna_group_info$n_cells_per_grna[["non-targeting"]]))

# 3. get covariate matrix
covariate_matrix_df <- response_odm |>
  get_cell_covariates() |>
  dplyr::slice(idxs) |>
  dplyr::select(n_nonzero, n_umis, bio_rep, phase, p_mito)
row.names(covariate_matrix_df) <- NULL
Z <- model.matrix(object = formula(~ log(n_nonzero) + log(n_umis) + bio_rep + phase + p_mito),
                                 data = covariate_matrix_df)
colnames(Z) <- c("intercept", "lg_n_nonzero", "lg_n_umis", "bio_rep_d1", "bio_rep_d2", "phase_d1", "phase_d2", "p_mito")

# 4. select gene
set.seed(3)
ex_gene <- response_odm |>
  get_feature_covariates() |>
  dplyr::arrange(desc(mean_expression)) |>
  dplyr::slice(1:100) |>
  dplyr::sample_n(1) |>
  row.names()
orig_y <- as.numeric(response_odm[[ex_gene, idxs]])

# 5. get the model for y | Z
fit_y_orig <- MASS::glm.nb(formula = orig_y ~ . + 0, data = as.data.frame(Z))
theta <- fit_y_orig$theta
mus_y <- as.numeric(fit_y_orig$fitted.values)

# 6. get the model for x | Z
fit_x <- glm(formula = orig_x ~ . + 0, family = binomial(), data = as.data.frame(Z))
mus_x <- as.numeric(fit_x$fitted.values)


#################################
# STEP 2: RUN SIMULATION FUNCTION
#################################
run_simulation <- function(Y, idx_mat, Z, theta_hypothesized, n_sim = NULL) {
  if (is.null(n_sim)) n_sim <- ncol(Y)
  sapply(seq(1, n_sim), function(i) {
    print(paste0("Running simulation ", i))
    y <- Y[,i]
    # regress the synthetic Y onto Z
    fit <- glm(y ~ Z + 0,
               family = MASS::negative.binomial(theta = theta_hypothesized))

    # extract z-scores
    z_scores <- sceptre2:::run_glm_perm_score_test_with_ingredients(Z = Z,
                                                                    working_resid = fit$residuals,
                                                                    w = fit$weights,
                                                                    index_mat = idx_mat)
    z_star <- z_scores[1]
    z_null <- z_scores[-1]
    z_star

    # compute theoretical and empirical p-values
    p_theory <- 2 * pnorm(q = -abs(z_star), lower.tail = TRUE)
    p_camp <- sceptre2:::compute_empirical_p_value(z_star = z_star, z_null = z_null, "both")

    # finally, compute permutation
    ts <- apply(X = idx_mat, MARGIN = 2, FUN = function(col) mean(y[col]))
    t_star <- ts[1]
    t_null <- ts[-1]
    p_perm <- sceptre2:::compute_empirical_p_value(t_star, t_null, "both")

    c(p_theory = p_theory, p_camp = p_camp, p_perm = p_perm)
  }) |> t()
}


##############################################
# STEP 3: GENERATE CORRELATED DATA AND RUN SIM
##############################################
n_sim <- 2000
# y first
Y <- sapply(X = mus_y, FUN = function(mu_y) MASS::rnegbin(n = n_sim, mu = mu_y, theta = theta)) |> t()

# generate x so that x is correlated with z
x <- sapply(X = mus_x, FUN = function(mu_x) rbinom(n = 1, size = 1, prob = mu_x))
x_idx <- which(x == 1)
B <- 100000
x_tilde <- replicate(n = B, expr = sample.int(n = length(x), size = sum(x))) - 1L
idx_mat <- cbind(matrix(x_idx, ncol = 1), x_tilde)

# run sim
sim_res_correlated <- run_simulation(Y = Y, idx_mat = idx_mat, Z = Z,
                                     theta_hypothesized = theta)

saveRDS(object = sim_res_correlated, file = paste0(result_dir, "correlated_sim_result.rds"))

########################################
# GENERATE UNCORRELATED DATA AND RUN SIM
########################################
# keep y from above
# generate x fresh
x <- rbinom(n = length(orig_x), size = 1, prob = mean(orig_x))
x_idx <- which(x == 1)
x_tilde <- replicate(n = B, expr = sample.int(n = length(x), size = sum(x))) - 1L
idx_mat <- cbind(matrix(x_idx, ncol = 1), x_tilde)

# run sim
sim_res_uncorrelated <- run_simulation(Y = Y, idx_mat = idx_mat, Z = Z, theta_hypothesized = 5 * theta, n_sim = 100)



#######
#
#######

as.data.frame(sim_res_correlated) |>
  tidyr::pivot_longer(cols = c("p_theory", "p_camp", "p_perm"),
                      names_to = "method", values_to = "p_value") |>
  ggplot(mapping = aes(y = p_value, col = method)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black")




######################
# ASSESS RUNNING TIME
######################
y <- Y[,1]
# 1. regress y on Z
glm_time <- microbenchmark(fit <- glm(y ~ Z + 0,
                                      family = MASS::negative.binomial(theta = theta)),
                           times = 30, unit = "s") |> summary()

# 2. permute x B times
x_idx <- which(x == 1)
B <- 100000
random_idx_time <- microbenchmark(x_tilde <- replicate(n = B, expr = sample.int(n = length(x), size = sum(x))),
                                  times = 5, unit = "s") |> summary()
idx_mat <- cbind(matrix(x_idx, ncol = 1), x_tilde) - 1L


# 3. get null z-scores
z_score_time <- microbenchmark(z_scores <- sceptre2:::run_glm_perm_score_test_with_ingredients(Z = Z,
                                                                                               working_resid = fit$residuals,
                                                                                               w = fit$weights,
                                                                                               index_mat = idx_mat), times = 30, unit = "s") |> summary()
z_star <- z_scores[1]
z_null <- z_scores[-1]
ks.test(z_scores, pnorm)
hist(z_scores)
abline(v = z_star)

p_theoretical <- 2 * pnorm(q = -abs(z_star), lower.tail = TRUE)
p_perm <- sceptre2:::compute_empirical_p_value(z_star = z_star, z_null = z_null, "both")

# put the times into a data frame
time_df <- data.frame(time = c(random_idx_time[["median"]], glm_time[["median"]], z_score_time[["median"]]),
                      operation = c("Generate permutation idxs", "Fit GLM", "Compute null statistics"))

