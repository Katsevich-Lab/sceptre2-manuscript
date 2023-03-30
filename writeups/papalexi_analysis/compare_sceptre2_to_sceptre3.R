LOCAL_SCEPTRE2_DATA_DIR <-.get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
papalexi_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "data/papalexi/eccite_screen/")

# gene info
gene_odm_fp <- paste0(papalexi_dir, "gene/matrix.odm")

# grna info
grna_odm_fp <- paste0(papalexi_dir, "grna_expression/matrix.odm")

# mm odm metadata fp
mm_metadata_fp <- paste0(papalexi_dir, "multimodal_metadata.rds")

# construct mm odm
mm_odm <- ondisc::read_multimodal_odm(odm_fps = c(gene_odm_fp, grna_odm_fp),
                                      multimodal_metadata_fp = mm_metadata_fp)

# get the in-memory gene matrix
gene_odm <- mm_odm |> ondisc::get_modality("gene")
response_matrix <- gene_odm[[seq(1, nrow(gene_odm)),]]
rownames(response_matrix) <- ondisc::get_feature_ids(gene_odm)

# get the in-memory grna matrix
grna_odm <- mm_odm |> ondisc::get_modality("grna_expression")
grna_matrix <- grna_odm[[seq(1, nrow(grna_odm)),]]
rownames(grna_matrix) <- ondisc::get_feature_ids(grna_odm)

# covariate matrix
covariate_data_frame <- mm_odm |> ondisc::get_cell_covariates()

# grna group data frame
grna_group_data_frame <- data.frame(grna_id = rownames(grna_odm@feature_covariates),
                                    grna_group = grna_odm@feature_covariates$target)

# set formulas, grna group target name
gene_formula <- ~ log(gene_n_umis) + log(gene_n_nonzero) + bio_rep + p_mito
grna_group <- "target"

# set hyperparameters
B <- 25000
side <- "both"
max_b_per_batch <- 25000
in_memory <- TRUE
statistic <- "full" # "distilled" is faster but might be less powerful
return_dist <- FALSE
screen_b <- 500

# select gene-gRNA group pairs to analyze
gene_grna_group_pairs <- expand.grid(response_id = mm_odm |>
                                       ondisc::get_modality("gene") |>
                                       ondisc::get_feature_ids(),
                                     grna_group = c("CUL3")) |> dplyr::sample_n(100)

# analyze the data using sceptre2
result_2 <- sceptre2::run_sceptre_low_moi(mm_odm = mm_odm,
                                          response_grna_group_pairs = gene_grna_group_pairs,
                                          form = gene_formula,
                                          response_modality_name = "gene",
                                          grna_modality_name = "grna_expression",
                                          grna_group_column_name = "target",
                                          B = B,
                                          side = side,
                                          max_b_per_batch = max_b_per_batch,
                                          in_memory = in_memory,
                                          statistic = statistic,
                                          return_dist = return_dist,
                                          screen_b = screen_b)


result_3 <- sceptre3::run_sceptre_lowmoi(response_matrix = response_matrix,
                                         grna_matrix = grna_matrix,
                                         covariate_data_frame = covariate_data_frame,
                                         grna_group_data_frame = grna_group_data_frame,
                                         formula_object = gene_formula,
                                         calibration_check = FALSE,
                                         response_grna_group_pairs = gene_grna_group_pairs)


res <- dplyr::left_join(x = result_3 |> dplyr::select(p_value, grna_group, response_id, round),
                        y = result_2, by = c("grna_group", "response_id"),
                        suffix = c("_sceptre_3", "_sceptre_2"))

library(ggplot2)
p1 <- ggplot(data = res, mapping = aes(x = p_value_sceptre_3,
                                       y = p_value_sceptre_2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw()

p2 <- ggplot(data = res, mapping = aes(x = p_value_sceptre_3,
                                       y = p_value_sceptre_2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  scale_x_continuous(trans = katlabutils::revlog_trans()) +
  scale_y_continuous(trans = katlabutils::revlog_trans())

ggsave(filename = "~/Desktop/sceptre_version_comparison_untrans.pdf", plot = p1, device = "pdf", scale = 0.9, width = 4, height = 4)
ggsave(filename = "~/Desktop/sceptre_version_comparison_trans.pdf", plot = p2, device = "pdf", scale = 0.9, width = 4, height = 4)
