library(ondisc)
sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- "gasperini"

FRAC_EXPRESSED_TRHESH <- 0.005
N_CELLS_PER_GRNA_THRESH <- 10

for (paper in papers) {
  paper_dir <- paste0(sceptre2_data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  for (dataset in datasets) {
    print(paste0("paper: ", paper, " dataset: ", dataset))
    mm_odm <- lowmoi::read_all_modalities(paper, dataset)

    # i. perform cell qc; restrict attention to cells with at least 1 grna and 1 gene expressed
    global_cell_covariates <- mm_odm |> get_cell_covariates()
    ok_cells <- global_cell_covariates |> dplyr::filter(grna_assignment_n_nonzero >= 1) |> row.names()
    mm_odm_sub <- mm_odm[,ok_cells]

    # ii. perform feature QC
    modalities <- names(mm_odm_sub@modalities)
    # grna assignment modality; keep features expressed in N_CELLS_PER_GRNA_THRESH cells
    grna_assign_modality <- get_modality(mm_odm_sub, "grna_assignment")
    n_cells_per_grna <- grna_assign_modality |>
      get_feature_covariates() |>
      dplyr::pull(n_nonzero)
    grna_assign_modality <- grna_assign_modality[n_cells_per_grna >= N_CELLS_PER_GRNA_THRESH,]
    mm_odm_sub@modalities$grna_assignment <- grna_assign_modality

    # grna expression modality: keep the same features as above
    grna_expression_modality <- get_modality(mm_odm_sub, "grna_expression")
    ok_grnas <- get_feature_ids(grna_assign_modality)
    grna_expression_modality <- grna_expression_modality[ok_grnas,]
    mm_odm_sub@modalities[["grna_expression"]] <- grna_expression_modality

    # response modalities; keep features expressed in FRAC_EXPRESSED_TRHESH of cells
    remaining_modalities <- modalities[!(modalities %in% c("grna_assignment", "grna_expression"))]
    for (modality in remaining_modalities) {
      modality_odm <- get_modality(mm_odm_sub, modality)
      feats_to_keep <- get_highly_expressed_features(modality_odm, FRAC_EXPRESSED_TRHESH)
      mm_odm_sub@modalities[[modality]] <- modality_odm[feats_to_keep,]
    }

    # iii. Write all modalities
    lowmoi::save_all_modalities(multimodal_odm = mm_odm_sub, paper = paper,
                                dataset = dataset, metadata_file_name = "metadata_qc.rds")

    # v. create a multimodal ondisc matrix free of redundancy and write
    mm_odm_sub_proc <- lowmoi::process_multimodal_odm(mm_odm_sub, FALSE)
    multimodal_metadata_fp <- paste0(paper_dir, dataset, "/multimodal_metadata.rds")
    save_multimodal_odm(multimodal_odm = mm_odm_sub_proc, multimodal_metadata_fp = multimodal_metadata_fp)
  }
}
