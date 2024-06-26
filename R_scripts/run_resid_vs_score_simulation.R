# simulation study # 2: comparing score statistic to mean over residuals statistic
# NB regression model with two covariates over draws, randomly vary (1) size parameter
# and (2) treatment coefficient (null vs. alternative).
args <- commandArgs(trailingOnly = TRUE)

library(camp)
library(tidyverse)
library(rlecuyer)
n_outer_rep <- 500
proc_id <- as.integer(args[1])

# seed setting
.lec.SetPackageSeed(4) |> invisible()
snames <- as.character(seq(1, n_outer_rep))
.lec.CreateStream(snames) |> invisible()
.lec.CurrentStream(snames[proc_id]) |> invisible()

# define the parameters that control the simulation
gamma <- c(-0.6, 0.8, 0.9)
n_rep <- 500
n_cells <- 5000
mu_z <- c(0.0, 0.0)
rho <- 0.5
frac_null <- 0.9

# generate the covariate matrix
z <- cbind(1, MASS::mvrnorm(n = n_cells,
                            mu = c(0.0, 0.0),
                            Sigma = toeplitz(c(1, rho))))
mu_x <- as.numeric(binomial()$linkinv(z %*% gamma))
x <- rbinom(n = n_cells, size = 1, prob = mu_x)
covariate_matrix <- cbind(z, x)
colnames(covariate_matrix) <- c("intercept", "z1", "z2", "x")
covariate_matrix_intercept_free <- covariate_matrix[,-1]
covariate_matrix_intercept_free_x_free <- covariate_matrix[,c(-1, -4)]

# generate permutation idxs
permutations <- permute_bernoulli_treatment_vector(x)

# generate results matrix
colnames <- c("p_resid", "p_score", "p_lrt", "null_true", "resid_time", "score_time", "lrt_time")
m <- as.data.frame(matrix(nrow = n_rep, ncol = 7))
colnames(m) <- colnames

# sample the runs under the alternative hypothesis
null_idxs <- sample(seq(1, n_rep), size = frac_null * n_rep) |> sort()

for (i in seq(1, n_rep)) {
  if (i %% 5 == 0) print(i)
  theta <- runif(1, 0.1, 5)
  null_true <- i %in% null_idxs
  beta <- c(0.9, 0.1, 0.3, if (null_true) 0.0 else 0.1)
  family_object <- MASS::negative.binomial(theta = theta)
  y <- generate_glm_data(
    design_matrix = covariate_matrix,
    coefficients = beta,
    family_object = family_object,
    add_intercept = FALSE
  )

  # fit the reduced GLM
  glm_fit_time <- system.time(
    fit_reduced <- glm(
      y ~ covariate_matrix_intercept_free_x_free,
      family = family_object,
    )
  )[["elapsed"]]

  # run residual precomputation and compute residual p-value
  resid_time <- system.time({
    precomputation_residual <- run_resid_precomputation(fit_reduced)
    p_resid <- run_perm_test_resid_stat_binary_trt(permutations, precomputation_residual)$p
  })[["elapsed"]]

  # run score precomputation and comptue score p-value
  score_time <- system.time({
    precomputation_score <- run_score_stat_precomputation(fit_reduced)
    p_score <- run_perm_test_score_stat_binary_trt(permutations, precomputation_score)$p
  })[["elapsed"]]

  # compute a standard GLM Wald p-value
  lrt_time <- system.time({
    fit_full <- glm(
      y ~ covariate_matrix_intercept_free,
      family = family_object,
    )
    lrt_test <- anova(fit_reduced, fit_full, test = "Chisq")
    p_lrt <- lrt_test$`Pr(>Chi)`[2]
  })[["elapsed"]]

  m[i,] <- c(p_resid = p_resid,
             p_score = p_score,
             p_lrt = p_lrt,
             null_true = null_true,
             resid_time = glm_fit_time + resid_time,
             score_time = glm_fit_time + score_time,
             lrt_time = glm_fit_time + lrt_time)
}

LOCAL_SCEPTRE2_DATA_DIR <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
dir_to_save <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "results/extra_analyses/resid_vs_score_sim/")
if (!dir.exists(dir_to_save)) dir.create(dir_to_save, recursive = TRUE)
saveRDS(object = m,
        file = paste0(dir_to_save, "raw_result_", proc_id, ".rds"))
