library(ondisc)
sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- list.files(sceptre2_data_dir)
papers <- papers[papers != "gasperini"]

# This script performs cell-wise QC, among other operations, on low MOI data.
# We (i) restrict attention to cells that received a single gRNA (as determined by the original authors) and
# (ii) filter for cells that passed other QC metrics implemented by the original authors (stored in the "passed_qc" column).

# set params
N_CELLS_PER_GRNA_THRESH <- 10
FRAC_EXPRESSED_TRHESH <- 0.005

# 0) save all modalities function
save_all_modalities <- function(multimodal_odm, paper, dataset, metadata_file_name, sceptre2_data_dir = paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")) {
  dataset_dir <- paste0(sceptre2_data_dir, paper, "/", dataset)
  modality_list <- multimodal_odm@modalities |> names()
  for (modality in modality_list) {
    save_odm(odm = get_modality(multimodal_odm, modality),
             metadata_fp = paste0(dataset_dir, "/", modality, "/", metadata_file_name))
  }
}

# 1) process multimodal ODM function
process_multimodal_odm <- function(mm_odm) {
  cell_covariate_m <- mm_odm |> get_cell_covariates()
  cell_covariate_m <- cell_covariate_m |> dplyr::mutate(grna_assignment_n_nonzero = NULL,
                                                        grna_assignment_n_umis = NULL)
  cell_covariate_colnames <- colnames(cell_covariate_m)
  shared_covariates <- c("batch", "p_mito", "phase", "bio_rep", "lane")
  for (shared_covariate in shared_covariates) {
    match_col_names <- grep(pattern = shared_covariate, x = cell_covariate_colnames, value = TRUE)
    if (length(match_col_names) >= 1) {
      cell_covariate_m[[shared_covariate]] <- cell_covariate_m[[match_col_names[1]]]
      for (match_col_name in match_col_names) cell_covariate_m[[match_col_name]] <- NULL
    }
  }
  mm_odm@global_cell_covariates <- cell_covariate_m
  return(mm_odm)
}

# 2) Set the MIMOSCA formula objects
mimosca_formula_objs <- list(frangieh = formula(~ n_nonzero + n_umis + phase + batch + 0),
                             schraivogel = formula(~ n_nonzero + n_umis + batch + 0),
                             papalexi = formula(~ n_nonzero + n_umis + bio_rep + phase + p_mito + 0),
                             liscovitch = formula(~ n_nonzero + n_fragments + 0),
                             simulated = formula(~ n_nonzero + n_umis + 0))

nb_regression_formula_objs <- list(frangieh = "~ offset(log(n_umis)) + log(n_nonzero) + phase + batch",
                                   schraivogel = "~ offset(log(n_umis)) + log(n_nonzero) + batch",
                                   papalexi = "~ offset(log(n_umis)) + log(n_nonzero) + bio_rep + phase + p_mito",
                                   liscovitch = "~ offset(log(n_fragments))",
                                   simulated = "~ offset(log(n_umis)) + log(n_nonzero)")

# 2) loop over datasets, loading all modalities
for (paper in papers) {
  paper_dir <- paste0(sceptre2_data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  for (dataset in datasets) {
    # load the dataset into a multimodal ODM
    print(paste0("paper: ", paper, " dataset: ", dataset))
    multimodal_metadata_fp <- paste0(paper_dir, dataset, "/multimodal_metadata.rds")
    if (file.exists(multimodal_metadata_fp)) file.remove(multimodal_metadata_fp)
    mm_odm <- lowmoi::read_all_modalities(paper, dataset)

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
    # grna assignment modality: keep features expressed in N_CELLS_PER_GRNA_THRESH cells. Also, add a "gRNA_assigned" column to the cell covariate matrix.
    grna_assign_modality <- get_modality(mm_odm_sub, "grna_assignment")
    grna_assign_mat <- lowmoi::load_whole_odm(grna_assign_modality)
    gRNA_assignments <- apply(X = grna_assign_mat,
                              MARGIN = 2,
                              FUN = function(col) names(which.max(col))) |> unname()
    grna_assign_modality <- grna_assign_modality |>
      mutate_cell_covariates(assigned_gRNA = gRNA_assignments)
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
      feats_to_keep <- get_highly_expressed_features(modality_odm, FRAC_EXPRESSED_TRHESH)
      mm_odm_sub@modalities[[modality]] <- modality_odm[feats_to_keep,]
    }

    # iii. perform feature ID cleanup; remove underscores and replace with dashes for all features
    for (modality in modalities) {
      modality_odm <- get_modality(mm_odm_sub, modality)
      modality_odm@ondisc_matrix@feature_ids <- gsub(pattern = "_", replacement = "-", x = modality_odm@ondisc_matrix@feature_ids, fixed = TRUE)
      row.names(modality_odm@feature_covariates) <- gsub(pattern = "_", replacement = "-", x = row.names(modality_odm@feature_covariates), fixed = TRUE)
      if (modality == "grna_assignment") {
        modality_odm <- mutate_cell_covariates(modality_odm, assigned_gRNA = gsub(pattern = "_", replacement = "-", x = assigned_gRNA, fixed = TRUE))
      }
      mm_odm_sub@modalities[[modality]] <- modality_odm
    }

    # iv. add the mimosca/nb formula objects to each response modality
    for (modality in remaining_modalities) {
      modality_odm <- get_modality(mm_odm_sub, modality)
      modality_odm@misc[["mimosca_formula"]] <- mimosca_formula_objs[[paper]]
      modality_odm@misc[["nb_regression_formula"]] <- nb_regression_formula_objs[[paper]]
      mm_odm_sub@modalities[[modality]] <- modality_odm
    }
    
    # Write the modified multimodal odm
    save_all_modalities(multimodal_odm = mm_odm_sub, paper = paper, dataset = dataset, metadata_file_name = "metadata_qc.rds")
    
    # v. create a multimodal ondisc matrix free of redundancy and write
    mm_odm_sub_proc <- process_multimodal_odm(mm_odm_sub)
    save_multimodal_odm(multimodal_odm = mm_odm_sub_proc, multimodal_metadata_fp = multimodal_metadata_fp)
  }
}
