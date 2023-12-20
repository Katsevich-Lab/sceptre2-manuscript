# simulation study # 2: comparing score statistic to mean over residuals statistic
# NB regression model with two covariates
# over draws, randomly vary (1) size parameter and (2) treatment coefficient (null vs. alternative)
set.seed(5)
library(camp)
library(katlabutils)
conflicts_prefer(dplyr::filter)
library(tidyverse)
library(cowplot)

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
m <- as.data.frame(matrix(nrow = n_rep, ncol = 6))
colnames(m) <- c("p_resid", "p_score", "p_lrt", "null_true", "resid_time", "score_time")

for (i in seq(1, n_rep)) {
  if (i %% 5 == 0) print(i)
  theta <- runif(1, 0.1, 5)
  null_true <- as.logical(rbinom(1, 1, frac_null))
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
  fit_full <- glm(
    y ~ covariate_matrix_intercept_free,
    family = family_object,
  )
  lrt_test <- anova(fit_reduced, fit_full, test = "Chisq")
  p_lrt <- lrt_test$`Pr(>Chi)`[2]
  m[i,] <- c(p_resid = p_resid, p_score = p_score, p_lrt = p_lrt,
             null_true = null_true, resid_time = glm_fit_time + resid_time,
             score_time = glm_fit_time + score_time)
}

# process result data frame and save
m <- m |>
  dplyr::mutate(p_resid_trans = -log(p_resid),
                p_score_trans = -log(p_score),
                p_lrt_trans = -log(p_lrt),
                null_true = (null_true == 1))

# apply bh
fdr_level <- 0.1
n_nonnull <- sum(!m$null_true)
summary_df <- m |> select(p_resid, p_score, p_lrt, null_true) |> 
  pivot_longer(cols = c("p_resid", "p_score", "p_lrt"),
               names_to = "method", values_to = "p_val") |>
  group_by(method) |>
  mutate(p_adj = p.adjust(p_val, method = "BH"), signif = p_adj < fdr_level) |>
  filter(signif) |>
  summarize(n_total_discoveries = dplyr::n(),
            n_false_discoveries = sum(null_true))

# assess mean running time
time_result_df <- m |> select(resid_time, score_time) |>
  pivot_longer(cols = c("resid_time", "score_time"),
               names_to = "method", values_to = "time") |>
  group_by(method) |>
  summarize(m_resid_time = mean(time), 
            m_score_time = mean(time),
            upper_ci = m_resid_time + 1.96 * sd(time)/sqrt(n_rep),
            lower_ci = m_resid_time - 1.96 * sd(time)/sqrt(n_rep))

f_p <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/extra_analyses/score_vs_resid_sim.rds")
to_save <- list(result_df = m, summary_df = summary_df, time_result_df = time_result_df)
saveRDS(object = to_save, file = f_p)