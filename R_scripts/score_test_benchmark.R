library(statmod)
library(camp)
library(dplyr)
library(tidyr)
set.seed(4)

# define variables describing experiment
prob_perts <- c(0.001, 0.01, 0.05, 0.1, 0.5)
n_rep_per_setting <- 50
n_trt_vect <- 1000
n <- 100000
family_object <- MASS::negative.binomial(5)

# sapply over prob_perts
res <- sapply(prob_perts, function(prob_pert) {
  # generate the treatment idxs
  m <- matrix(data = rbinom(n_trt_vect * n, size = 1, prob_pert),
              nrow = n, ncol = n_trt_vect)
  trt_idxs <- apply(X = m, MARGIN = 2, FUN = function(col) which(col == 1) - 1L)
  
  # iterate over reps per setting
  out <- sapply(seq(1, n_rep_per_setting), FUN = function(rep_id) {
    z <- MASS::mvrnorm(n = n, mu = c(-0.5, 0.5), Sigma = toeplitz(c(1, 0.5)))
    y <- generate_glm_data(
      design_matrix = z,
      coefficients = c(0.6, 0.1, 0.3),
      family_object = family_object,
      add_intercept = TRUE
    )
    # fit the model
    fit <- glm(y ~ z, family = family_object)
    
    # compute the z-scores using eigen decomp
    eigen_time <- system.time({
      precomputation_score <- run_score_stat_precomputation(fit)
      for (i in seq(1, n_trt_vect)) {
        z_eigen <- camp:::compute_observed_full_statistic(a = precomputation_score$a,
                                                          w = precomputation_score$w,
                                                          D = precomputation_score$D,
                                                          s = length(trt_idxs[[i]]),
                                                          trt_idxs = trt_idxs[[i]])   
      }
    })
    
    # compute the z-scores using statmod
    statmod_time <- system.time({
      statmod_z <- statmod::glm.scoretest(fit = fit, x2 = m, dispersion = 1)
    })
      
    # output results
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
