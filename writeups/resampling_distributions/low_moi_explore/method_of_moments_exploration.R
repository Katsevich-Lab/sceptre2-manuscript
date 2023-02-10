library(evmix)

samp_norm <- rnorm(250000, sd = 1)
samp <- samp_norm[samp_norm > 2] - 2
n <- length(samp)

# gpd mle
fit_mle <- fgpd(x = samp)
xi_mle <- fit_mle$mle[2]
sigmau_mle <- fit_mle$mle[1]

# gpd mm
m <- mean(samp)
v <- var(samp)
xi_mm <- (1/2) * (1 - m^2/v)
sigmau_mm <- m * (1 - xi_mm)

# pareto mle
samp_p1 <- samp + 1
alpha_hat <- length(samp_p1)/(sum(log(samp_p1)))


xgrid <- seq(0, max(samp), length.out = 501)
y_mle <- dgpd(xgrid, sigmau = sigmau_mle, xi = xi_mle)
y_mm <- dgpd(xgrid, sigmau = sigmau_mm, xi = xi_mm)

hist(samp, freq = FALSE, breaks = 15, ylim = c(0, max(y_mle, y_mm)))
lines(x = xgrid, y = y_mle, col = "blue")
lines(x = xgrid, y = y_mm, col = "red")

cdf_funct_factory <- function(xi, sigma) {
  f <- function(x) {
    ifelse(test = xi < 0 & x > -sigma/xi, 1e-16, 1 - (1 + xi * x/sigma)^(-1/xi) )
  }
  return(f)
}

cdf <- cdf_funct_factory(xi_mle, sigmau_mle)
curve(cdf, from = 0, to = 4)
