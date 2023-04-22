#############
# SCHRAIVOGEL
#############
response_odm <- load_dataset_modality("schraivogel/enhancer_screen_chr11/gene")
grna_odm <- load_dataset_modality("schraivogel/enhancer_screen_chr11/grna_assignment")

# strategy 1: get target assignments via max op
grna_targets <- get_target_assignments_via_max_op(grna_odm)
trt_vect <- response_odm[["ZFPM2", grna_targets == "ZFPM2-enh"]] |> as.numeric()
cntrl_vect <- response_odm[["ZFPM2", grna_targets == "non-targeting"]] |> as.numeric()

sum(trt_vect >= 1) # fails to pass QC!
sum(cntrl_vect >= 1)

# strategy 2: get target assignments via sceptre
grna_matrix <- grna_odm[[seq(1, nrow(grna_odm)),]]
feature_df <- grna_odm |> ondisc::get_feature_covariates()
rownames(grna_matrix) <- rownames(feature_df)
grna_group_df <- data.frame(grna_id = rownames(feature_df),
                            grna_group = feature_df$target)
grna_assignments <- sceptre:::assign_grnas_to_cells_lowmoi_v2(grna_matrix = grna_matrix,
                                          grna_group_data_frame = grna_group_df,
                                          calibration_check = FALSE,
                                          n_calibration_pairs = NULL)
trt_vect <- response_odm[["ZFPM2", grna_assignments$grna_group_idxs$`ZFPM2-enh`]]
cntrl_vect <- response_odm[["ZFPM2", grna_assignments$all_nt_idxs]]
sum(trt_vect >= 1)
sum(cntrl_vect >= 1)

# strategy 3: call the sceptre method
response_mat <- response_odm[[seq(1, nrow(response_odm)),]]
rownames(response_mat) <- ondisc::get_feature_ids(response_odm)
covariate_data_frame <- response_odm |> ondisc::get_cell_covariates()

res <- sceptre::run_sceptre_lowmoi(response_matrix = response_mat,
                            grna_matrix = grna_matrix,
                            covariate_data_frame = covariate_data_frame,
                            grna_group_data_frame = grna_group_df,
                            formula_object = formula(~n_nonzero + n_umis),
                            response_grna_group_pairs = data.frame(response_id = c("ZFPM2", "AKIP1"),
                                                                   grna_group = c("ZFPM2-enh", "ZFPM2-enh")),
                            calibration_check = FALSE)

# ALL AGREE: SOMETHING IS WRONG WITH THE SAMPLE SIZE DF OR THE FUNCTION TO COMPUTE SAMPLE SIZES USING THIS DF
response_odm <- load_dataset_modality("schraivogel/enhancer_screen_chr11/gene")
grna_odm <- load_dataset_modality("schraivogel/enhancer_screen_chr11/grna_assignment")

##########
# PAPALEXI
##########
response_odm <- load_dataset_modality("papalexi/eccite_screen/gene")
grna_odm <- load_dataset_modality("papalexi/eccite_screen/grna_assignment")
expression_vector <- response_odm[["WIPF3",]] |> as.numeric()
response_mat <- response_odm[[seq(1, nrow(response_odm))]]
rownames(response_mat) <- ondisc::get_feature_ids(response_odm)
expression_vector_2 <- response_mat["WIPF3",]
response_mat <- sceptre:::set_matrix_accessibility(matrix_in = response_mat, make_row_accessible = TRUE)
expression_vector_3 <- sceptre:::load_csr_row(j = response_mat@j,
                                              p = response_mat@p,
                                              x = response_mat@x,
                                              row_idx = which(rownames(response_mat) == "WIPF3"),
                                              n_cells = ncol(response_mat))

# strategy 2: get target assignments via sceptre
#grna_matrix <- grna_odm[[seq(1, nrow(grna_odm)),]]
#feature_df <- grna_odm |> ondisc::get_feature_covariates()
#rownames(grna_matrix) <- rownames(feature_df)
#grna_group_df <- data.frame(grna_id = rownames(feature_df),
#                            grna_group = feature_df$target)
#grna_assignments <- sceptre:::assign_grnas_to_cells_lowmoi_v2(grna_matrix = grna_matrix,
#                                                              grna_group_data_frame = grna_group_df,
#                                                              calibration_check = TRUE,
#                                                              n_calibration_pairs = NULL)
#nt1_idxs <- grna_assignments$all_nt_idxs[grna_assignments$indiv_nt_grna_idxs$NTg1]
