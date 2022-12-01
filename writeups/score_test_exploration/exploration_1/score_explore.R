##########
# FULL GLM
##########
run_experiment <- function(fitted_means = FALSE, n_rep = 2000, correlated = FALSE) {
  n <- 5000
  beta_v <- c(1, 2, 1)
  dat <- matrix(data = c(rep(1, n),
                         x1 = rbinom(n, 1, 0.4),
                         x2 = rnorm(n)),
                ncol = 3)
  lin_pred <- as.numeric(dat %*% beta_v)
  mus <- exp(lin_pred)
  
  if (correlated) {
    x_test <- sapply(binomial()$linkinv(dat[,3] - 0.1), function(curr_mu) rbinom(n = 1, size = 1, prob = curr_mu))
  } else {
    x_test <- rbinom(n = n, size = 1, prob = 0.5)
  }
  
  out <- sapply(X = seq(1, n_rep), FUN = function(i) {
    print(i)
    y <- sapply(mus, function(curr_mu) rpois(n = 1, lambda = curr_mu))
    
    # we do a good job of estimating the parameters via glm
    fit <- glm(y ~ dat[,-1], family = poisson)
    mu_hat <- fit$fitted.values
    z_score_full <- statmod::glm.scoretest(fit = fit, x_test)
    
    # next, get the distilled z-score
    y_1 <- y[x_test == 1]
    mu_hat_1 <- if (fitted_means) mu_hat[x_test == 1] else mus[x_test == 1]
    s_y <- sum(y_1)
    s_mu_hat <- sum(mu_hat_1)
    
    # get the distilled z-score
    z_score_distilled <- (s_y - s_mu_hat)/sqrt(s_mu_hat)
    
    # get the normalized pearson residual sum
    pearson_resid <- (y_1 - mu_hat_1)/sqrt(mu_hat_1)
    sum_pearson <- sum(pearson_resid)
    norm_sum_pearson <- (1/sqrt(length(pearson_resid))) * sum_pearson

    # get the studentized residuals
    #h <- VGAM::hatvalues(fit)[x_test == 1]
    #boot_resid <- pearson_resid/sqrt(1 - h)
    #norm_boot_resid <- (1/sqrt(length(boot_resid))) * sum(boot_resid)
    
    # get the gcm
    r <- y_1 - mu_hat_1
    gcm <- 1/sqrt(length(r)) * sum(r)/sd(r)
    
    c(z_score_distilled = z_score_distilled,
      z_score_full = z_score_full,
      norm_sum_pearson = norm_sum_pearson,
      gcm = gcm)
    
    # finally, get the empirically-normalized residual
    # emp_norm_raw_resid <- (1/sqrt(n)) * sum(r_1)/sd(r_1)
    # return
    # c(z_score_distilled = z_score_distilled,
    #  norm_sum_pearson = norm_sum_pearson,
    #  emp_norm_raw_resid = emp_norm_raw_resid,
    #  z_score_full = z_score_full)
  }) |> t()  
}

##################
# Pt 1: true means
##################

# true means
out_true_means <- run_experiment(fitted_means = FALSE, n_rep = 2000)
x_grid <- seq(from = -3.5, to = 3.5, by = 0.01)
y <- dnorm(x_grid)

hist(out_true_means[,"z_score_full"], freq = FALSE, breaks = 20, main = "Full z-score")
lines(x = x_grid, y = y, col = "red")

hist(out_true_means[,"z_score_distilled"], freq = FALSE, breaks = 20, main = "Distilled z-score")
lines(x = x_grid, y = y, col = "red")

hist(out_true_means[,"norm_sum_pearson"], freq = FALSE,
     breaks = 20, main = "Normalized sum of Pearson residuals")
lines(x = x_grid, y = y, col = "red")

hist(out_true_means[,"gcm"], freq = FALSE,
     breaks = 20, main = "GCM statistic")
lines(x = x_grid, y = y, col = "red")

####################
# Pt 2: fitted means
####################
out_fitted_means <- run_experiment(fitted_means = TRUE, n_rep = 2000, correlated = FALSE)

hist(out_fitted_means[,"z_score_full"], freq = FALSE,
     breaks = 20, main = "Full z-score")
lines(x = x_grid, y = y, col = "red")

hist(out_fitted_means[,"z_score_distilled"], freq = FALSE,
     breaks = 20, main = "Distilled z-score")
lines(x = x_grid, y = y, col = "red")

hist(out_fitted_means[,"norm_sum_pearson"], freq = FALSE,
     breaks = 20, main = "Normalized sum of Pearson residuals")
lines(x = x_grid, y = y, col = "red")

hist(out_fitted_means[,"gcm"], freq = FALSE,
     breaks = 20, main = "gcm")
lines(x = x_grid, y = y, col = "red")

plot(out_fitted_means[,"z_score_full"], out_fitted_means[,"z_score_distilled"])
abline(a = 0, b = 1, col = "red")

####################################
# Pt 3 fitted means with correlation
####################################
out_fitted_means <- run_experiment(fitted_means = TRUE, n_rep = 2000, correlated = TRUE)

hist(out_fitted_means[,"z_score_full"], freq = FALSE,
     breaks = 20, main = "Full z-score")
lines(x = x_grid, y = y, col = "red")

hist(out_fitted_means[,"z_score_distilled"], freq = FALSE,
     breaks = 20, main = "Distilled z-score")
lines(x = x_grid, y = y, col = "red")

hist(out_fitted_means[,"norm_sum_pearson"], freq = FALSE,
     breaks = 20, main = "Normalized sum of Pearson residuals")
lines(x = x_grid, y = y, col = "red")

hist(out_fitted_means[,"gcm"], freq = FALSE,
     breaks = 20, main = "gcm")
lines(x = x_grid, y = y, col = "red")

plot(out_fitted_means[,"z_score_full"], out_fitted_means[,"z_score_distilled"])
abline(a = 0, b = 1, col = "red")

###################
# Pt 3: comparisons
###################
plot(out_fitted_means[,"z_score_distilled"], out_fitted_means[,"z_score_full"])
abline(a = 0, b = 1, col = "red")
plot(out_fitted_means[,"z_score_distilled"], out_fitted_means[,"norm_sum_pearson"])
abline(a = 0, b = 1, col = "red")

# The distilled z-score seems to differ from the full z-score by a constant factor.
# Meanwhile, the distilled 
# At any rate, these three quantities are NOT the same. (And we of course have used the canonical link here.)
