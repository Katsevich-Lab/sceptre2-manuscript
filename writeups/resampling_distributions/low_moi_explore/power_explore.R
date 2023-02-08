library(eva)
library(tidyverse)
library(katlabutils)
library(cowplot)

sceptre2_results_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
# load empirical distributions
resampling_dists <- readRDS(paste0(sceptre2_results_dir, "resampling_distributions/sceptre_resampling_dists.rds"))


# gpd_mat <- list()
# samples_list <- c(seq(1e4, 1e3, length.out = 10))
# for (k in 1:length(samples_list)) {
#   gpd_mat[[k]] <- matrix(0, ncol = 330, nrow = 10)
#   for (i in 1:330) {
#     if (i %in% c(163, 177, 198)){
#       next
#     }
#     working_dir <- getwd()
#     subDir <- sprintf("figures/power_exploration/%d_resamples_%.2f_percentile",  samples_list[k], 0.96)
#     gpd_param <- read.csv(sprintf("%s/%d.csv", file.path(working_dir, subDir), i))
#     gpd_mat[[k]][1:5, i] <- gpd_param[, i+1]
#     gpd_mat[[k]][6:10, i] <- gpd_param[, i+331]
#   }
#   print(k)
# }

# k-s test using fitted distribution
samples_list <- c(round(seq(5e4, 1e3, length.out = 10)))
power_mat <- list()
for (i in 1:length(samples_list)){
  power_mat[[i]] <- matrix(0, nrow = 1000, ncol = 330)
  no_resamp <- samples_list[i]
  q <- 0.96
  for (r in 73:100) {
    for (l in 1:330){
      # set.seed
      set.seed(r)
      idx <- l
      samples <- sample(1:5e5, no_resamp)
      z_null <- as.numeric(resampling_dists[idx, 4+samples])
      
      for(tail in c("left", "right")){
        sign_flip <- if(tail == "right") 1 else -1
        z_null_pos <- sign_flip*z_null
        u_pos <- quantile(z_null_pos, q)
        u <- sign_flip*u_pos
        
        # extract the data
        trun_data <- z_null_pos[z_null_pos > u_pos]
        n <- sum(z_null_pos > u_pos)
        
        # method of moment
        xi_mom <- (1 - (mean(trun_data) - u_pos)^2/var(trun_data))/2
        sig_mom <- (mean(trun_data) - u_pos)*(1 - xi_mom)
        
        sim_data <- rgpd(1e6, loc = u_pos, scale = sig_mom, shape = xi_mom)
        
        # k-s test
        test_ks <- ks.test(trun_data, sim_data)
        
        # store p_value and test statistic
        if(tail == "right"){
          power_mat[[i]][6, r] <- u_pos
          power_mat[[i]][7, r] <- sig_mom
          power_mat[[i]][8, r] <- xi_mom
          power_mat[[i]][9,r] <- test_ks$statistic
          power_mat[[i]][10,r] <- test_ks$p.value
        }else{
          power_mat[[i]][1, r] <- u_pos
          power_mat[[i]][2, r] <- sig_mom
          power_mat[[i]][3, r] <- xi_mom
          power_mat[[i]][4,r] <- test_ks$statistic
          power_mat[[i]][5,r] <- test_ks$p.value
        }
      }   
      print(power_mat[[i]][1:4,r])
    }
    print(r)
  }
}

