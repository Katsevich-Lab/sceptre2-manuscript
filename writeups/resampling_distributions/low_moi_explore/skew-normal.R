library(tidyverse)
library(sn)

sceptre2_results_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
# load empirical distributions
resampling_dists <- readRDS(paste0(sceptre2_results_dir, "resampling_distributions/sceptre_resampling_dists.rds"))

idx <- 1
z_null <- as.numeric(resampling_dists[idx,5:ncol(resampling_dists)])
tibble(z_null) |>
  ggplot(aes(x = z_null)) + 
  geom_histogram(aes(y = ..density..), bins = 1000) + 
  stat_function(fun = function(x)(dsn(x, dp = out@param$dp))) + 
  coord_cartesian(xlim = c(1, 3))

out <- sn::selm(z_null ~ 1, family = "SN")

ks.test(z_null, "psn", dp = out@param$dp)


dsn(x = 3, dp = out@param$dp)

out@dp