library(ondisc) # devtools::install_github('timothy-barry/ondisc')
#devtools::install_github('timothy-barry/sceptre3')
library(sceptre3)
library(BH)

LOCAL_SCEPTRE2_DATA_DIR <-.get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
papalexi_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "data/papalexi/eccite_screen/")

# gene info
gene_odm_fp <- paste0(papalexi_dir, "gene/matrix.odm")
gene_metadata_fp <- paste0(papalexi_dir, "gene/metadata_qc.rds")
gene_odm <- read_odm(odm_fp = gene_odm_fp, metadata_fp = gene_metadata_fp)
gene_covariate_matrix <- gene_odm |> get_cell_covariates() 
gene_expression_matrix <- gene_odm[[seq(1, nrow(gene_odm)),]]
rownames(gene_expression_matrix) <- get_feature_ids(gene_odm)

# grna info
grna_odm_fp <- paste0(papalexi_dir, "grna_expression/matrix.odm")
grna_metadata_fp <- paste0(papalexi_dir, "grna_expression/metadata_qc.rds")
grna_odm <- read_odm(odm_fp = grna_odm_fp, metadata_fp = grna_metadata_fp)
grna_matrix <- grna_odm[[seq(1, nrow(grna_odm)),]]
grna_groups <- data.frame(grna_id = rownames(grna_odm@feature_covariates),
                          grna_group = grna_odm@feature_covariates$target)

# protein info
protein_odm_fp <- paste0(papalexi_dir, "protein/matrix.odm")
protein_metadata_fp <- paste0(papalexi_dir, "protein/metadata_qc.rds")
protein_odm <- read_odm(odm_fp = protein_odm_fp, metadata_fp = protein_metadata_fp)
protein_covariate_matrix <- protein_odm |> get_cell_covariates()
protein_matrix <- protein_odm[[seq(1,nrow(protein_odm)),]]
rownames(protein_matrix) <- get_feature_ids(protein_odm)

# set formulas, grna group target name
gene_formula <- ~ log(n_umis) + log(n_nonzero) + bio_rep + p_mito
protein_formula <- ~ log(n_umis) + bio_rep + p_mito

#######################################
# SET ARGS FOR GENE EXPRESSION ANALYSIS
#######################################
response_matrix <- gene_expression_matrix
grna_matrix <- grna_matrix
rownames(grna_matrix) <- ondisc::get_feature_ids(grna_odm)
covariate_data_frame <- gene_covariate_matrix
grna_group_data_frame <- grna_groups
formula_object <- gene_formula
calibration_check <- FALSE
unique_grna = unique(grna_groups$grna_group)
response_grna_group_pairs <- expand.grid(response_id = get_feature_ids(gene_odm),
                                           grna_group = unique_grna[-which(unique_grna == 'non-targeting')]) # an example set of pairs

test_stat <- "full"
return_resampling_dist <- FALSE
adaptive_permutation_test <- TRUE
fit_skew_normal <- TRUE

result_gene <- run_sceptre_lowmoi(response_matrix = response_matrix,
                                  grna_matrix = grna_matrix,
                                  covariate_data_frame = covariate_data_frame,
                                  grna_group_data_frame = grna_group_data_frame,
                                  formula_object = formula_object,
                                  response_grna_group_pairs = response_grna_group_pairs,
                                  calibration_check = calibration_check)

##########################################
# SET ARGS FOR PROTEIN EXPRESSION ANALYSIS
##########################################
response_matrix <- as.matrix(protein_matrix)
grna_matrix <- grna_matrix
rownames(grna_matrix) <- ondisc::get_feature_ids(grna_odm)
covariate_data_frame <- protein_covariate_matrix
grna_group_data_frame <- grna_groups
formula_object <- protein_formula
calibration_check <- FALSE
response_grna_group_pairs <- expand.grid(response_id = rownames(response_matrix),
                                         grna_group = unique_grna[-which(unique_grna == 'non-targeting')]) # an example set of pairs
test_stat <- "full"
return_resampling_dist <- FALSE
adaptive_permutation_test <- TRUE
fit_skew_normal <- TRUE

result_protein <- run_sceptre_lowmoi(response_matrix = response_matrix,
                                     grna_matrix = grna_matrix,
                                     covariate_data_frame = covariate_data_frame,
                                     grna_group_data_frame = grna_group_data_frame,
                                     formula_object = formula_object,
                                     response_grna_group_pairs = response_grna_group_pairs,
                                     calibration_check = calibration_check)

sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
output_filename <- paste0(sceptre2_dir, "results/papalexi_analysis/")
saveRDS(result_gene,paste0(output_filename,"sceptre_full_mrna_results_with_effect_size.rds"))
saveRDS(result_protein,paste0(output_filename,"sceptre_full_protein_results_with_effect_size.rds"))

