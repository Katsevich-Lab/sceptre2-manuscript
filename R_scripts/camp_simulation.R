library(ondisc)
library(sceptre2)
library(microbenchmark)
library(ggplot2)
library(katlabutils)
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
ex_gene <- "CXCL10"
orig_y <- as.numeric(response_odm[[ex_gene, idxs]])

# 5. get the model for y | Z
fit_y_orig <- MASS::glm.nb(formula = orig_y ~ . + 0, data = as.data.frame(Z))
y_coef <- coef(fit_y_orig)
theta <- fit_y_orig$theta
mus_y <- exp(as.numeric(Z %*% y_coef))

# 6. get the model for x | Z
b <- binomial()
fit_x <- glm(formula = orig_x ~ . + 0, family = binomial(), data = as.data.frame(Z))
x_coef <- coef(fit_x)
mus_x <- b$linkinv(as.numeric(Z %*% x_coef))

#################################
# STEP 2: RUN SIMULATION FUNCTION
#################################
run_simulation <- function(Y, idx_mat, Z, theta_hypothesized, n_sim = NULL, return_null_dist = FALSE, approx = TRUE) {
  resamp_dist <- list()
  if (is.null(n_sim)) n_sim <- ncol(Y)
  
  out_m <- matrix(nrow = n_sim, ncol = 3)
    for (i in seq(1, n_sim)) {
      print(paste0("Running simulation ", i))
      y <- Y[,i]
      
      # regress the synthetic Y onto Z
      fit <- glm(y ~ Z + 0,
                 family = MASS::negative.binomial(theta = theta_hypothesized))
      
      # extract z-scores
      if (i == 1 || !approx) {
        z_scores <- sceptre2:::run_glm_perm_score_test_with_ingredients(Z = Z,
                                                                        working_resid = fit$residuals,
                                                                        w = fit$weights,
                                                                        index_mat = idx_mat)
        z_star <- z_scores[1]
        z_null <- z_scores[-1]
      } else {
        z_star <- sceptre2:::run_glm_perm_score_test_with_ingredients(Z = Z,
                                                                      working_resid = fit$residuals,
                                                                      w = fit$weights,
                                                                      index_mat = idx_mat[,1,drop = FALSE])
      }

      # compute theoretical and empirical p-values
      p_theory <- 2 * pnorm(q = -abs(z_star), lower.tail = TRUE)
      p_camp <- sceptre2:::compute_empirical_p_value(z_star = z_star, z_null = z_null, "both")
      
      # finally, compute permutation
      ts <- sceptre2:::low_level_permutation_test(y = y, index_mat = idx_mat)
      t_star <- ts[1]
      t_null <- ts[-1]
      p_perm <- sceptre2:::compute_empirical_p_value(t_star, t_null, "both")
      out_m[i,] <- c(p_theory = p_theory, p_camp = p_camp, p_perm = p_perm)
      if (i == 1) {
        resamp_dist[["camp_null"]] <- z_null
        resamp_dist[["camp_star"]] <- z_star
        resamp_dist[["perm_null"]] <- ts
        resamp_dist[["perm_star"]] <- t_star
      }
    }
    colnames(out_m) <- c("p_theory", "p_camp", "p_perm")
    return(list(out_m = out_m, resamp_dist = resamp_dist))
}

##############################################
# STEP 3: GENERATE CORRELATED DATA AND RUN SIM
##############################################
# n_sim <- 2000
n_sim <- 5
set.seed(3)
# y first
Y <- sapply(X = mus_y, FUN = function(mu_y) MASS::rnegbin(n = n_sim, mu = mu_y, theta = theta)) |> t()

# keep x fixed
x_idx <- which(orig_x == 1)
B <- 100000
x_tilde <- replicate(n = B, expr = sample.int(n = length(orig_x), size = sum(orig_x)))
idx_mat <- cbind(matrix(x_idx, ncol = 1), x_tilde) - 1L

# run sim
sim_res_correlated_corret_model <- run_simulation(Y = Y, idx_mat = idx_mat,
                                                  Z = Z, theta_hypothesized = theta,
                                                  return_null_dist = TRUE, approx = FALSE)

###################################################
# STEP 4: CORRELATED AND MODEL MISSPECIFICATION SIM
###################################################
sim_res_correlated_misspec <- run_simulation(Y = Y, idx_mat = idx_mat, Z = Z,
                                             theta_hypothesized =  5 * theta,
                                             return_null_dist = TRUE, approx = FALSE)

#########################################
# STEP 6: UNCORRELATED, GLM INCORRECT SIM
#########################################
# keep y from above
# generate x fresh
x <- rbinom(n = length(orig_x), size = 1, prob = mean(orig_x))
x_idx <- which(x == 1)
x_tilde <- replicate(n = B, expr = sample.int(n = length(x), size = sum(x)))
idx_mat <- cbind(matrix(x_idx, ncol = 1), x_tilde) - 1L

# run sim
sim_res_uncorrelated_misspec <- run_simulation(Y = Y, idx_mat = idx_mat, Z = Z,
                                       theta_hypothesized = 5 * theta,
                                       return_null_dist = TRUE, approx = FALSE)

#######################################
# STEP 7: UNCORRELATED, GLM CORRECT SIM
#######################################
# run sim
sim_res_uncorrelated_correct_model <- run_simulation(Y = Y, idx_mat = idx_mat, Z = Z,
                                                     theta_hypothesized = theta,
                                                     return_null_dist = TRUE, approx = FALSE)

# save
sim_res_all <- list(sim_res_correlated_corret_model = sim_res_correlated_corret_model,
     sim_res_correlated_misspec = sim_res_correlated_misspec,
     sim_res_uncorrelated_correct_model = sim_res_uncorrelated_correct_model,
     sim_res_uncorrelated_misspec = sim_res_uncorrelated_misspec)

saveRDS(object = sim_res_all,
        file = paste0(result_dir, "simulation_study_res.rds"))

######################
# ASSESS RUNNING TIME
######################
if (FALSE) {
  y <- Y[,1]
  # 1. regress y on Z
  glm_time <- microbenchmark(fit <- glm(y ~ Z + 0,
                                        family = MASS::negative.binomial(theta = theta)),
                             times = 30, unit = "s") |> summary()
  
  # 2. permute x B times
  x_idx <- which(x == 1)
  B <- 25000
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
}
