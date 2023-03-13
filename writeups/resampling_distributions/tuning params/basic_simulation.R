xi <- 0.5
omega <- 1
alpha <- 3
n <- 4999
y <- c(sn::rsn(n = n, xi = xi, omega = omega, alpha = alpha),
       sn::rsn(n = n, xi = 0.5 * xi, omega = -3 * omega, alpha = 2 * alpha))

fit_and_evaluate_skew_normal(-5, y)
compute_empirical_p_value(null_statistics = y, z_orig = -5, side = 0L)

hist(y, freq = FALSE, breaks = 100)
fit <- fit_skew_normal_funct(y)
x_grid <- seq(from = range(y)[1], to = range(y)[2], length.out = 100)
y_grid <- sn::dsn(x = x_grid, xi = fit[1], omega = fit[2], alpha = fit[3])
lines(x_grid, y_grid, col = "blue")
abline(v = quantile(y, .9), col = "red")
abline(v = quantile(y, .995), col = "darkgreen")

ret <- check_sn_tail(y, fit[1], fit[2], fit[3])
ret
