##########
# FULL GLM
##########
n <- 5000
beta_v <- c(1, 2, 1)
dat <- matrix(data = c(rep(1, n),
              x1 = rbinom(n, 1, 0.4),
              x2 = rnorm(n)),
              ncol = 3)
lin_pred <- as.numeric(dat %*% beta_v)
mus <- exp(lin_pred)
x_test <- rbinom(n = n, size = 1, prob = 0.5)

out <- sapply(X = seq(1, 2000), FUN = function(i) {
  print(i)
  y <- sapply(mus, function(curr_mu) rpois(n = 1, lambda = curr_mu))
  
  # we do a good job of estimating the parameters via glm
  fit <- glm(y ~ dat[,-1], family = poisson)
  mu_hat <- fit$fitted.values
  z_score_full <- statmod::glm.scoretest(fit = fit, x_test)
  
  # next, get the distilled z-score
  y_1 <- y[x_test == 1]
  exp_o_1 <- fit$fitted.values[x_test == 1]
  # exp_o_1 <- mus[x_test == 1]
  s_y <- sum(y_1)
  s_exp_o <- sum(exp_o_1)
  
  # get the distilled z-score
  z_score_distilled <- (s_y - s_exp_o)/sqrt(s_exp_o)
  
  c(z_score_distilled = z_score_distilled,
    z_score_full = z_score_full)
  
  # get the normalized pearson residual sum
  # r_1 <- y_1 - exp_o_1
  # sum_pearson <- sum((y_1 - exp_o_1)/sqrt(exp_o_1))
  # norm_sum_pearson <- (1/sqrt(n)) * sum_pearson
  # finally, get the empirically-normalized residual
  # emp_norm_raw_resid <- (1/sqrt(n)) * sum(r_1)/sd(r_1)
  
  # return
  # c(z_score_distilled = z_score_distilled,
  #  norm_sum_pearson = norm_sum_pearson,
  #  emp_norm_raw_resid = emp_norm_raw_resid,
  #  z_score_full = z_score_full)
}) |> t()

# There clearly appears to be a difference between 

z_score_distilled <- out[,"z_score_distilled"]
z_score_full <- out[,"z_score_full"]
lm(z_score_full ~ z_score_distilled)

plot(z_score_distilled, z_score_full)
abline(0, 1, col = "red")

x_grid <- seq(from = -3.5, to = 3.5, by = 0.01)
y <- dnorm(x_grid)
hist(z_score_full, freq = FALSE, ylim = c(0, 0.4), breaks = 20)
lines(x = x_grid, y = y, col = "red")

x_grid <- seq(from = -3.5, to = 3.5, by = 0.01)
y <- dnorm(x_grid, 0, 0.7)
hist(z_score_distilled, freq = FALSE, ylim = c(0, 0.6), breaks = 20)
lines(x = x_grid, y = y, col = "red")

mean(z_score_full)
sd(z_score_full)

mean(z_score_distilled)
sd(z_score_distilled)

############################
# SIMPLE MEAN + OFFSET MODEL
############################
o <- runif(n, 2, 4)
exp_o <- exp(o)

z <- replicate(1000, {
  y <- sapply(exp_o, function(curr_exp_o) rpois(n = 1, lambda = curr_exp_o))
  s_y <- sum(y)
  s_exp_o <- sum(exp_o)
  z_score_distilled <- (s_y - s_exp_o)/sqrt(s_exp_o)
  z_score_distilled
})

hist(z, freq = FALSE, ylim = c(0, 0.4))
lines(x = x_grid, y = y, col = "red")
