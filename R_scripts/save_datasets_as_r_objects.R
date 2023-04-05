#####################################################################################
# This script saves the portions of the Papalexi and Frangieh IFN-gamma datasets
# that are used by the run_sceptre_lowmoi function for the discovery analysis (Fig 5)
# in standard R format (as opposed to ondisc format). This script thus facillitates
# replication of the full discovery analyses.
#####################################################################################

####################
# FRANGIEH IFN-GAMMA
####################
LOCAL_SCEPTRE2_DATA_DIR <-.get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
ifn_gamma_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "data/frangieh/ifn_gamma/")

# gene info
gene_odm_fp <- paste0(ifn_gamma_dir, "gene/matrix.odm")

# grna info
grna_odm_fp <- paste0(ifn_gamma_dir, "grna_assignment/matrix.odm")

# mm odm metadata fp
mm_metadata_fp <- paste0(ifn_gamma_dir, "multimodal_metadata.rds")

# construct mm odm
mm_odm <- ondisc::read_multimodal_odm(odm_fps = c(gene_odm_fp, grna_odm_fp),
                                      multimodal_metadata_fp = mm_metadata_fp)

# get the in-memory gene matrix
gene_odm <- mm_odm |> ondisc::get_modality("gene")
response_matrix <- gene_odm[[seq(1, nrow(gene_odm)),]]
rownames(response_matrix) <- ondisc::get_feature_ids(gene_odm)

# get the in-memory grna matrix
grna_odm <- mm_odm |> ondisc::get_modality("grna_assignment")
grna_matrix <- grna_odm[[seq(1, nrow(grna_odm)),]]
rownames(grna_matrix) <- ondisc::get_feature_ids(grna_odm)

# covariate matrix
covariate_data_frame <- mm_odm |> ondisc::get_cell_covariates() |>
  dplyr::select(gene_n_nonzero, gene_n_umis)

# grna group data frame
grna_group_data_frame <- data.frame(grna_id = rownames(grna_odm@feature_covariates),
                                    grna_group = grna_odm@feature_covariates$target)

# set formulas, grna group target name
formula_object <- mm_odm@global_misc$formula

# set the gene-grna group pairs
response_grna_group_pairs <- sceptre::generate_all_pairs(response_matrix = response_matrix,
                                                         grna_group_data_frame = grna_group_data_frame)

# create the list of items to write
l <- list(response_matrix = response_matrix,
          grna_matrix = grna_matrix,
          covariate_data_frame = covariate_data_frame,
          grna_group_data_frame = grna_group_data_frame,
          formula_object = formula_object, 
          response_grna_group_pairs = response_grna_group_pairs)

dir_to_save <- paste0(ifn_gamma_dir)
file_to_save <- paste0(dir_to_save, "r_objects.rds")
saveRDS(object = l, file = file_to_save)

###############
# PAPALEXI GENE
###############
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
covariate_data_frame <- mm_odm |> ondisc::get_cell_covariates() |>
  dplyr::select(gene_n_umis, gene_n_nonzero, bio_rep, p_mito)

# grna group data frame
grna_group_data_frame <- data.frame(grna_id = rownames(grna_odm@feature_covariates),
                                    grna_group = grna_odm@feature_covariates$target)

# set formulas, grna group target name
formula_object <- mm_odm@global_misc$formula

# set the gene-grna group pairs
response_grna_group_pairs <- sceptre::generate_all_pairs(response_matrix = response_matrix,
                                                         grna_group_data_frame = grna_group_data_frame)

# create list
l <- list(response_matrix = response_matrix,
          grna_matrix = grna_matrix,
          covariate_data_frame = covariate_data_frame,
          grna_group_data_frame = grna_group_data_frame,
          formula_object = formula_object,
          response_grna_group_pairs = response_grna_group_pairs)
dir_to_save <- paste0(papalexi_dir)
file_to_save <- paste0(dir_to_save, "r_objects.rds")
saveRDS(object = l, file = file_to_save)
