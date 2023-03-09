#########################################
# Read and process ground truth results
#########################################

undershoot <- read_csv("undershoot_res_fit.csv")[,-1]
overshoot <- read_csv("overshoot_res_fit.csv")[,-1]
quantile_list <- seq(0.03, 0.97, length.out = 10)
no_sam <- round(exp(seq(log(1e3), log(5e4), length.out = 10)))
# rearrange the data frame
B <- 100

undershoot_df_right <- data.frame(id = rep(1:B, 10*10),
                                  ratio_value = 0,
                                  no_sam = 0,
                                  ratio_quantile = 0)
overshoot_df_right <- data.frame(id = rep(1:B, 10*10),
                                 ratio_value = 0,
                                 no_sam = 0,
                                 ratio_quantile = 0)
undershoot_df_left <- data.frame(id = rep(1:B, 10*10),
                                 ratio_value = 0,
                                 no_sam = 0,
                                 ratio_quantile = 0)
overshoot_df_left <- data.frame(id = rep(1:B, 10*10),
                                ratio_value = 0,
                                no_sam = 0,
                                ratio_quantile = 0)


# i: quantile; j: no of sample
for (i in 1:10) {
  for (j in 1:10) {
    start <- (j - 1 + (i-1)*10)*B +1
    end <- (j + (i-1)*10)*B
    undershoot_df_right[start:end, 2] <- as.vector(undershoot[(((j-1)*B+1) :(j*B)), (i-1)*3+32])[[1]]
    undershoot_df_right[start:end, 3] <- rep(no_sam[j], B)
    undershoot_df_right[start:end, 4] <- rep(quantile_list[i], B)
    overshoot_df_right[start:end, 2] <- as.vector(overshoot[(((j-1)*B+1) :(j*B)), (i-1)*3+32])[[1]]
    overshoot_df_right[start:end, 3] <- rep(no_sam[j], B)
    overshoot_df_right[start:end, 4] <- rep(quantile_list[i], B)
    
    undershoot_df_left[start:end, 2] <- as.vector(undershoot[(((j-1)*B+1) :(j*B)), (i-1)*3+2])[[1]]
    undershoot_df_left[start:end, 3] <- rep(no_sam[j], B)
    undershoot_df_left[start:end, 4] <- rep(quantile_list[i], B)
    overshoot_df_left[start:end, 2] <- as.vector(overshoot[(((j-1)*B+1) :(j*B)), (i-1)*3+2])[[1]]
    overshoot_df_left[start:end, 3] <- rep(no_sam[j], B)
    overshoot_df_left[start:end, 4] <- rep(quantile_list[i], B)
  }
}

subsampling_results_ground_truth <- bind_rows(
  undershoot_df_left |> mutate(tail = "left", error_type = "undershoot"),
  undershoot_df_right |> mutate(tail = "right", error_type = "undershoot"),
  overshoot_df_left |> mutate(tail = "left", error_type = "overshoot"),
  overshoot_df_right |> mutate(tail = "right", error_type = "overshoot")
) |>
  mutate(evaluation = "ground truth") |>
  rename(ratio = ratio_value, num_samples_fit = no_sam) |>
  as_tibble()

#########################################
# Read and process estimated results
#########################################

undershoot <- read_csv("undershoot_res_power.csv")[,-1]
overshoot <- read_csv("overshoot_res_power.csv")[,-1]
quantile_list <- seq(0.03, 0.97, length.out = 10)
no_sam <- round(exp(seq(log(1e3), log(5e4), length.out = 10)))
# rearrange the data frame
B <- 100

undershoot_df_right <- data.frame(id = rep(1:B, 10*10),
                                  ratio_value = 0,
                                  no_sam = 0,
                                  ratio_quantile = 0)
overshoot_df_right <- data.frame(id = rep(1:B, 10*10),
                                 ratio_value = 0,
                                 no_sam = 0,
                                 ratio_quantile = 0)
undershoot_df_left <- data.frame(id = rep(1:B, 10*10),
                                 ratio_value = 0,
                                 no_sam = 0,
                                 ratio_quantile = 0)
overshoot_df_left <- data.frame(id = rep(1:B, 10*10),
                                ratio_value = 0,
                                no_sam = 0,
                                ratio_quantile = 0)


# i: quantile; j: no of sample
for (i in 1:10) {
  for (j in 1:10) {
    start <- (j - 1 + (i-1)*10)*B +1
    end <- (j + (i-1)*10)*B
    undershoot_df_right[start:end, 2] <- as.vector(undershoot[(((j-1)*B+1) :(j*B)), (i-1)*3+32])[[1]]
    undershoot_df_right[start:end, 3] <- rep(no_sam[j], B)
    undershoot_df_right[start:end, 4] <- rep(quantile_list[i], B)
    overshoot_df_right[start:end, 2] <- as.vector(overshoot[(((j-1)*B+1) :(j*B)), (i-1)*3+32])[[1]]
    overshoot_df_right[start:end, 3] <- rep(no_sam[j], B)
    overshoot_df_right[start:end, 4] <- rep(quantile_list[i], B)
    
    undershoot_df_left[start:end, 2] <- as.vector(undershoot[(((j-1)*B+1) :(j*B)), (i-1)*3+2])[[1]]
    undershoot_df_left[start:end, 3] <- rep(no_sam[j], B)
    undershoot_df_left[start:end, 4] <- rep(quantile_list[i], B)
    overshoot_df_left[start:end, 2] <- as.vector(overshoot[(((j-1)*B+1) :(j*B)), (i-1)*3+2])[[1]]
    overshoot_df_left[start:end, 3] <- rep(no_sam[j], B)
    overshoot_df_left[start:end, 4] <- rep(quantile_list[i], B)
  }
}

subsampling_results_estimate <- bind_rows(
  undershoot_df_left |> mutate(tail = "left", error_type = "undershoot"),
  undershoot_df_right |> mutate(tail = "right", error_type = "undershoot"),
  overshoot_df_left |> mutate(tail = "left", error_type = "overshoot"),
  overshoot_df_right |> mutate(tail = "right", error_type = "overshoot")
) |>
  mutate(evaluation = "estimate") |>
  rename(ratio = ratio_value, num_samples_fit = no_sam) |>
  as_tibble()

subsampling_results <- bind_rows(
  subsampling_results_ground_truth,
  subsampling_results_estimate
) |>
  rename(replicate = id) |>
  mutate(id = as.numeric(as.factor(ratio_quantile))) |>
  select(-ratio_quantile)

saveRDS(subsampling_results, "subsampling_nc_results.rds")
