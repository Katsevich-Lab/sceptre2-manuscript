library(statmod)
library(camp)
library(dplyr)
library(tidyr)
set.seed(4)

# define a set of sparsity levels
prob_perts <- c(0.001, 0.01, 0.05, 0.1, 0.5)
n_rep_per_setting <- 50

# sapply over prob_perts
res <- sapply(prob_perts, function(prob_pert) {
  n_rep <- 1000
  n <- 100000
  x <- rbinom(n = n, size = 1, prob = prob_pert)
  z <- MASS::mvrnorm(n = n, mu = c(-0.5, 0.5), Sigma = toeplitz(c(1, 0.5)))
  family_object <- MASS::negative.binomial(5)
  design_matrix <- cbind(x, z)
  y <- generate_glm_data(
    design_matrix = design_matrix,
    coefficients = c(0.6, 0.1, 0.1, 0.3),
    family_object = family_object,
    add_intercept = TRUE
  )
  # fit the model
  fit <- glm(
    y ~ design_matrix[,-x],
    family = family_object,
  )
  # obtain the permuted treatment vectors
  trt_idxs <- which(x == 1) - 1L
  s <- length(trt_idxs)
  m <- matrix(data = rep(x, n_rep), ncol = 1000)
  # compute the z-score using eigen decomp
  out <- sapply(seq(1, n_rep_per_setting), FUN = function(rep_id) {
    eigen_time <- system.time({
      precomputation_score <- run_score_stat_precomputation(fit)
      for (i in seq(1, n_rep)) {
        z_eigen <- camp:::compute_observed_full_statistic(a = precomputation_score$a,
                                                          w = precomputation_score$w,
                                                          D = precomputation_score$D,
                                                          s = s,
                                                          trt_idxs = trt_idxs)   
      }
    })
    # compute the z-score using statmod (i.e., QR decomp)
    statmod_time <- system.time({
      statmod_z <- statmod::glm.scoretest(fit = fit, x2 = m, dispersion = 1)
    })
    c(statmod_time = statmod_time[["elapsed"]],
      eigen_time = eigen_time[["elapsed"]],
      prob_pert = prob_pert,
      rep_id = rep_id)
  }) |> t() |> as.data.frame()
}, simplify = FALSE)
res_df <- res |> data.table::rbindlist()

# for each level, compute the mean runtime, as well as an upper and lower ci
res_df_sum <- res_df |>
  pivot_longer(cols = c("statmod_time", "eigen_time"), names_to = "method", values_to = "time") |>
  group_by(prob_pert, method) |>
  summarize(m_time = mean(time),
            lower_ci = m_time - 1.96/sqrt(n_rep_per_setting) * sd(time),
            upper_ci = m_time + 1.96/sqrt(n_rep_per_setting) * sd(time)) |>
  mutate(speedup = m_time[method == "statmod_time"]/m_time[method == "eigen_time"],
         sparsity = 100 * (1 - prob_pert))
res_df_sum
