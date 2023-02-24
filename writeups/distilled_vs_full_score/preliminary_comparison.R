data(gene_matrix_lowmoi)
data(grna_matrix_lowmoi)
data(covariate_data_frame_lowmoi)
data(grna_group_data_frame_lowmoi)
data(response_grna_group_pairs_lowmoi)

response_matrix <- gene_matrix_lowmoi
grna_matrix <- grna_matrix_lowmoi
covariate_data_frame <- covariate_data_frame_lowmoi
grna_group_data_frame <- grna_group_data_frame_lowmoi
calibration_check <- FALSE
formula_object <- formula(~log(gene_n_umis) + log(gene_n_nonzero) + bio_rep + phase + p_mito)
response_grna_group_pairs <- response_grna_group_pairs_lowmoi
test_stat <- "distilled"
return_resampling_dist <- FALSE
adaptive_permutation_test <- TRUE
fit_skew_normal <- TRUE

result_dist <- run_sceptre_lowmoi(response_matrix,
                                  grna_matrix,
                                  covariate_data_frame,
                                  grna_group_data_frame,
                                  formula_object,
                                  calibration_check,
                                  response_grna_group_pairs,
                                  test_stat,
                                  return_resampling_dist,
                                  adaptive_permutation_test,
                                  fit_skew_normal)

result_full <- run_sceptre_lowmoi(response_matrix,
                                  grna_matrix,
                                  covariate_data_frame,
                                  grna_group_data_frame,
                                  formula_object,
                                  calibration_check,
                                  response_grna_group_pairs,
                                  "full",
                                  return_resampling_dist,
                                  adaptive_permutation_test,
                                  fit_skew_normal,
                                  B2 = 25001L)

x <- dplyr::left_join(x = result_full, y = result_dist, by = c("response_id", "grna_group"), suffix = c("_full", "_dist"))

ggplot(data = x, mapping = aes(x = p_value_full, y = p_value_dist)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, col = "red")


p_dist_adj <- p.adjust(p = x$p_value_dist, method = "BH")
p_full_adj <- p.adjust(p = x$p_value_full, method = "BH")

# BH adjustment at FDR 0.1
sum(p_dist_adj < 0.1)
sum(p_full_adj < 0.1)

# BH adjustment at FDR 0.01
sum(p_dist_adj < 0.01)
sum(p_full_adj < 0.01)

ggplot(data = data.frame(p_dist_adj = p_dist_adj, p_full_adj = p_full_adj),
       mapping = aes(x = p_full_adj, y = p_dist_adj)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_vline(xintercept = 0.01, col = "blue") +
  geom_hline(yintercept = 0.01, col = "blue") +
  xlab("Full p-value (adjusted)") +
  ylab("Distilled p-valued (adjusted)")
