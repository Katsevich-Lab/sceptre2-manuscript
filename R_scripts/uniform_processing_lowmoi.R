library(ondisc)
sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- c("frangieh", "liscovitch",  "papalexi", "schraivogel", "simulated")

# This script performs cell-wise QC, among other operations, on low MOI data.
# We (i) restrict attention to cells that received a single grna (as determined by the original authors) and
# (ii) filter for cells that passed other QC metrics implemented by the original authors (stored in the "passed_qc" column).

# set params
FRAC_EXPRESSED_TRHESH <- 0.005

# 1.i) Set the MIMOSCA formula objects
mimosca_formula_objs <- list(frangieh = formula(~ n_nonzero + n_umis + phase + batch + 0),
                             schraivogel = formula(~ n_nonzero + n_umis + batch + 0),
                             papalexi = formula(~ n_nonzero + n_umis + bio_rep + phase + p_mito + 0),
                             liscovitch = formula(~ n_nonzero + n_fragments + 0),
                             simulated = formula(~ n_nonzero + n_umis + 0))

mimosca_formula_objs_protein <- list(frangieh = formula(~ n_umis + phase + batch + 0),
                                     papalexi = formula(~ n_umis + bio_rep + phase + p_mito + 0))


# 1.ii) Set the NB regression formula objects
nb_regression_formula_objs <- list(frangieh = "~ offset(log(n_umis)) + log(n_nonzero) + phase + batch",
                                   schraivogel = "~ offset(log(n_umis)) + log(n_nonzero) + batch",
                                   papalexi = "~ offset(log(n_umis)) + log(n_nonzero) + bio_rep + phase + p_mito",
                                   liscovitch = "~ offset(log(n_fragments))",
                                   simulated = "~ offset(log(n_umis)) + log(n_nonzero)")

nb_regression_formula_objs_protein <- list(frangieh = "~ offset(log(n_umis)) + phase + batch",
                                           papalexi = "~ offset(log(n_umis)) + bio_rep + phase + p_mito")

# 1.iii) Set the sceptre formula objects
sceptre_formula_objs <- list(frangieh = ~ log(response_n_umis) + log(response_n_nonzero) + phase + batch,
                             schraivogel = ~ log(response_n_umis) + log(response_n_nonzero) + batch,
                             papalexi = ~ log(response_n_umis) + log(response_n_nonzero) + bio_rep + phase + p_mito,
                             liscovitch = ~ log(response_n_fragments),
                             simulated = ~ log(response_n_umis) + log(response_n_nonzero))

sceptre_formula_objs_protein <- list(frangieh = ~ log(response_n_umis) + phase + batch,
                                     papalexi = ~ log(response_n_umis) + bio_rep + phase + p_mito)

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

    # i. perform cell QC; restrict attention to 1 grna/cell and "passed_qc" cells (if applicable)
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
    # grna assignment modality: keep features expressed in N_CELLS_PER_GRNA_THRESH cells. Also, add a "grna_assigned" column to the cell covariate matrix.
    grna_assign_modality <- get_modality(mm_odm_sub, "grna_assignment")
    grna_assign_mat <- lowmoi::load_whole_odm(grna_assign_modality)
    assigned_grna <- apply(X = grna_assign_mat,
                              MARGIN = 2,
                              FUN = function(col) names(which.max(col))) |> unname()
    grna_assign_modality <- grna_assign_modality |>
      mutate_cell_covariates(assigned_grna = assigned_grna)
    if (paper == "schraivogel") {
      grna_assign_modality <- grna_assign_modality |>
        mutate_feature_covariates(target = ifelse(is.na(known_effect), target, known_effect),
                                  known_effect = NULL)
    }

    # grna expression modality (if applicable): keep the same features as above
    if ("grna_expression" %in% modalities) {
      grna_expression_modality <- get_modality(mm_odm_sub, "grna_expression")
      if (paper == "schraivogel") {
        grna_expression_modality <- grna_expression_modality |>
          mutate_feature_covariates(target = known_effect, known_effect = NULL)
      }
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
        modality_odm <- mutate_cell_covariates(modality_odm, assigned_grna = gsub(pattern = "_", replacement = "-", x = assigned_grna, fixed = TRUE))
      }
      mm_odm_sub@modalities[[modality]] <- modality_odm
    }

    # iv. add the mimosca/nb formula objects to each response modality
    for (modality in remaining_modalities) {
      modality_odm <- get_modality(mm_odm_sub, modality)
      if (modality == "protein") {
        modality_odm@misc[["mimosca_formula"]] <- mimosca_formula_objs_protein[[paper]]
        modality_odm@misc[["nb_regression_formula"]] <- nb_regression_formula_objs_protein[[paper]]
        modality_odm@misc[["sceptre_formula"]] <- sceptre_formula_objs_protein[[paper]] 
      } else {
        if (dataset %in% c("enhancer_screen_chr11", "enhancer_screen_chr8") && paper == "schraivogel") { # special case: dataset is schraivogel/enhancer_screen_chr11 or schraivogel/enhancer_screen_chr8
          modality_odm@misc[["mimosca_formula"]] <- formula(~ n_nonzero + n_umis + 0)
          modality_odm@misc[["nb_regression_formula"]] <- "~ offset(log(n_umis)) + log(n_nonzero)"
          modality_odm@misc[["sceptre_formula"]] <- formula(~ log(response_n_umis) + log(response_n_nonzero))
        } else {
          modality_odm@misc[["mimosca_formula"]] <- mimosca_formula_objs[[paper]]
          modality_odm@misc[["nb_regression_formula"]] <- nb_regression_formula_objs[[paper]]
          modality_odm@misc[["sceptre_formula"]] <- sceptre_formula_objs[[paper]]  
        }
      }
      mm_odm_sub@modalities[[modality]] <- modality_odm
    }

    # Write all modalities
    lowmoi::save_all_modalities(multimodal_odm = mm_odm_sub, paper = paper, dataset = dataset, metadata_file_name = "metadata_qc.rds")

    # v. create a multimodal ondisc matrix free of redundancy and write
    mm_odm_sub_proc <- lowmoi::process_multimodal_odm(mm_odm_sub)
    save_multimodal_odm(multimodal_odm = mm_odm_sub_proc,
                        multimodal_metadata_fp = multimodal_metadata_fp)
    
    # vi. write the positive control pairs
    if (paper %in% c("frangieh", "papalexi", "schraivogel")) {
      # grouped pairs
      grna_assignment_modality <- mm_odm_sub_proc |> get_modality("grna_assignment")
      gene_modality <- mm_odm_sub_proc |> get_modality("gene")
      
      grna_feature_df <- grna_assignment_modality |> ondisc::get_feature_covariates()
      targets <- intersect(grna_feature_df |> dplyr::pull(target),
                           gene_modality |> ondisc::get_feature_ids())
      pc_pairs <- data.frame(grna_group = targets, response_id = targets)
      saveRDS(pc_pairs, file = paste0(paper_dir, dataset, "/pos_control_pairs_grouped.rds"))
      
      # ungrouped pairs
      ungroup_map <- data.frame(grna_id = row.names(grna_feature_df),
                                grna_group = grna_feature_df$target)
      ungroup_pc_pairs <- dplyr::left_join(ungroup_map, pc_pairs, by = "grna_group") |>
        na.omit() |> dplyr::select(grna_id, response_id)
      saveRDS(ungroup_pc_pairs, file = paste0(paper_dir, dataset, "/pos_control_pairs_single.rds"))
    }
  }
}
