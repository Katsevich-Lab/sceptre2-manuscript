library(dplyr)

n <- 5000
B <- 250000
theta <- 15

x1 <- rbinom(n, 1, 0.05)
x2 <- rnorm(n)

beta_v <- c(1, -0.6, 0.4, -0.3, 0.3)
Z <- matrix(data = c(intercept = rep(1, n),
                     x1 = rbinom(n, 1, 0.4),
                     x2 = rnorm(n),
                     x3 = x1 + x2 + rnorm(n),
                     x4 = x1 - 2 * x2 - rnorm(n)),
            ncol = 5)
lin_pred <- as.numeric(Z %*% beta_v)
mus <- exp(lin_pred)
y <- sapply(mus, function(curr_mu) MASS::rnegbin(n = 1, mu = curr_mu, theta = theta))
fit <- glm(y ~ Z + 0, family = MASS::neg.bin(theta))
index_mat <- replicate(n = B,
                       expr = sample.int(n = length(y), size = sum(x1))) - 1L
system.time(z_scores <- run_glm_perm_score_test(fit, index_mat)) # 0.5 s

X <- apply(X = index_mat[,1:1000], MARGIN = 2, FUN = function(col) {
  out <- integer(n)
  out[col + 1] <- 1
  return(out)
})

system.time(z_scores_wald <- apply(X, MARGIN = 2, function(col) {
  fit <- glm(y ~ Z + col + 0, family = MASS::neg.bin(theta))
  s <- summary(fit)
  s$coefficients["col", "t value"]
}))
(23.852 * 250)/(60^2) # 1 h 40 min
