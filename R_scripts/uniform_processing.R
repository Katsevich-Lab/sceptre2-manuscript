library(ondisc)
sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- list.files(sceptre2_data_dir)

# This script performs cell-wise QC.
# We (i) restrict attention to cells that received a single gRNA (as determined by the original authors) and
# (ii) filter for cells that passed other QC metrics implemented by the original authors (stored in the "passed_qc" column).

# set params
N_CELLS_PER_GRNA_THRESH <- 10
FRAC_EXPRESSED_TRHESH <- 0.005

# 0) General save and read multimodal_odm functions
# save
save_multimodal_odm <- function(multimodal_odm, paper, dataset, metadata_file_name, sceptre2_data_dir = paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")) {
  dataset_dir <- paste0(sceptre2_data_dir, paper, "/", dataset)
  modality_list <- multimodal_odm@modalities |> names()
  for (modality in modality_list) {
    save_odm(odm = get_modality(multimodal_odm, modality),
             metadata_fp = paste0(dataset_dir, "/", modality, "/", metadata_file_name))
  }
}

# read
read_multimodal_odm <- function(paper, dataset, sceptre2_data_dir = paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")) {
  dataset_dir <- paste0(sceptre2_data_dir, paper, "/", dataset)
  modality_vect <- list.files(dataset_dir)
  odm_list <- list()
  for (modality in modality_vect) {
    odm_dir <- paste0(dataset_dir, "/", modality, "/")
    curr_odm <- read_odm(odm_fp = paste0(odm_dir, "matrix.odm"),
                         metadata_fp = paste0(odm_dir, "metadata_orig.rds"))
    odm_list <- c(odm_list, curr_odm)
  }
  names(odm_list) <- modality_vect
  ret <- multimodal_ondisc_matrix(odm_list)
  return(ret)
}

# 1) Set the MIMOSCA formula objects
mimosca_formula_objs <- list(frangieh = formula(~ n_umis + phase + 0),
                             schraivogel = formula(~ n_umis + batch + 0),
                             papalexi = formula(~ n_umis + batch + phase + p_mito + 0),
                             liscovitch = formula(~ n_fragments + 0),
                             simulated = formula(~ n_umis + 0))

# 2) loop over datasets, loading all modalities
for (paper in papers) {
  paper_dir <- paste0(sceptre2_data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  for (dataset in datasets) {
    # load the dataset into a multimodal ODM
    print(paste0("paper: ", paper, " dataset: ", dataset))
    mm_odm <- read_multimodal_odm(paper, dataset)

    # i. perform cell QC; restrict attention to 1 gRNA/cell and "passed_qc" cells (if applicable)
    global_cell_covariates <- mm_odm |> get_cell_covariates()
    cell_logical_v <- global_cell_covariates$grna_assignment_n_nonzero == 1
    passed_qc_v <- grepl(pattern = "passed_qc", x = colnames(global_cell_covariates))
    if (any(passed_qc_v)) {
      passed_qc <- global_cell_covariates[, which(passed_qc_v)[1]]
      cell_logical_v <- cell_logical_v & passed_qc
    }
    mm_odm_sub <- mm_odm[,cell_logical_v]

    # ii. perform feature QC
    modalities <- names(mm_odm_sub@modalities)
    # grna assignment modality: keep features expressed in N_CELLS_PER_GRNA_THRESH cells
    grna_assign_modality <- get_modality(mm_odm_sub, "grna_assignment")
    grna_assign_mat <- lowmoi::load_whole_odm(grna_assign_modality)
    n_cells_per_gRNA <- Matrix::rowSums(grna_assign_mat)
    grnas_to_keep <- n_cells_per_gRNA >= N_CELLS_PER_GRNA_THRESH
    mm_odm_sub@modalities[["grna_assignment"]] <- grna_assign_modality[grnas_to_keep,]

    # grna expression modality (if applicable): keep the same features as above
    if ("grna_expression" %in% modalities) {
      grna_expression_modality <- get_modality(mm_odm_sub, "grna_expression")
      mm_odm_sub@modalities[["grna_expression"]] <- grna_expression_modality[grnas_to_keep,]
    }

    # response modalities: keep features expressed in FRAC_EXPRESSED_TRHESH of cells
    remaining_modalities <- modalities[!(modalities %in% c("grna_assignment", "grna_expression"))]
    for (modality in remaining_modalities) {
      modality_odm <- get_modality(mm_odm_sub, modality)
      exp_mat <- lowmoi::load_whole_odm(modality_odm)
      frac_expressed <- Matrix::rowSums(exp_mat >= 1)/ncol(exp_mat)
      feats_to_keep <- frac_expressed > FRAC_EXPRESSED_TRHESH
      mm_odm_sub@modalities[[modality]] <- modality_odm[feats_to_keep,]
    }

    # iii. perform feature ID cleanup; remove underscores and replace with dashes for all features
    for (modality in modalities) {
      modality_odm <- get_modality(mm_odm_sub, modality)
      modality_odm@ondisc_matrix@feature_ids <- gsub(pattern = "_", replacement = "-", x = modality_odm@ondisc_matrix@feature_ids, fixed = TRUE)
      row.names(modality_odm@feature_covariates) <- gsub(pattern = "_", replacement = "-", x = row.names(modality_odm@feature_covariates), fixed = TRUE)
      mm_odm_sub@modalities[[modality]] <- modality_odm
    }

    # iv. add the mimosca formula object to each response modality
    for (modality in remaining_modalities) {
      modality_odm <- get_modality(mm_odm_sub, modality)
      modality_odm@misc[["mimosca_formula"]] <- mimosca_formula_objs[[paper]]
      mm_odm_sub@modalities[[modality]] <- modality_odm
    }
    # Finally, write the modified multimodal odm
    save_multimodal_odm(multimodal_odm = mm_odm_sub, paper = paper, dataset = dataset, metadata_file_name = "metadata_qc.rds")
  }
}
