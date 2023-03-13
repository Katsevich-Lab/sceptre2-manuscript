mle_param_nc <- read_csv("figures/power_exploration/sknorm_tail_prob_500000_resamples_0.96_percentile/param_twosides.csv")
mom_param_nc <- read_csv("figures/power_exploration/sknorm_tail_prob_500000_resamples_0.96_percentile/mom_param.csv")

mle_param_twosides <- t(mle_param_nc[,-1])
mom_param_twosides <- t(mom_param_nc[,-1])
mle_overshoot_ratio <- as.numeric(mle_param_twosides[, 6])
mle_undershoot_ratio <- as.numeric(mle_param_twosides[, 7])
mom_overshoot_ratio <- as.numeric(mom_param_twosides[, 6])
mom_undershoot_ratio <- as.numeric(mom_param_twosides[, 7])

# accuracy using 5e5 mom estimate
mle_groundtruth <- data.frame(ratio = c(mle_overshoot_ratio, mle_undershoot_ratio), type = c(rep("overshoot", 660), rep("undershoot", 660)),
                              tail = c(rep(c(rep("left", 330),rep("right", 330)) , 2)), id = c(rep(1:330, 2)))
mom_groundtruth <- data.frame(ratio = c(mom_overshoot_ratio, mom_undershoot_ratio), type = c(rep("overshoot", 660), rep("undershoot", 660)),
                              tail = c(rep(c(rep("left", 330),rep("right", 330)) , 2)), id = c(rep(1:330, 2)))

groundtruth_results <- bind_rows(
  mle_groundtruth |> mutate(method = "MLE"),
  mom_groundtruth |> mutate(method = "MoM")
) |> 
  rename(error_type = type) |>
  as_tibble()

saveRDS(groundtruth_results, "ground_truth_nc_results.rds")
