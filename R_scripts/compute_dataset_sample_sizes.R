# for each dataset, comptue the following:
# i) for each gRNA, the number of cells receiving that gRNA, as well as the number of cells with nonzero gene expression
# i) furthermore, the type of each gRNA
# Thus, we seek to create a data frame with the following columns:
# i) dataset,
# ii) gRNA,
# iii) gene
# iii) n cells with gRNA
# iv) n with gRNA with nonzero gene expression
# v) gRNA type

library(ondisc)
sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
sceptre2_sample_sizes_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/dataset_sample_sizes/")
if (!dir.exists(sceptre2_sample_sizes_dir)) dir.create(sceptre2_sample_sizes_dir)
papers <- c("frangieh", "liscovitch",  "papalexi", "schraivogel", "simulated")

#######################################################
# FIRST: compute number of nonzero cells cells per gene
#######################################################
df <- lapply(papers, function(paper) {
  print(paste0("paper: ", paper))
  paper_dir <- paste0(sceptre2_data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  lapply(X = datasets, FUN = function(dataset) {
    print(paste0("paper: ", paper, " dataset: ", dataset))
    paper_fp <- paste0(paper, "/", dataset)
    mm_odm <- lowmoi::load_dataset_multimodal(paper_fp = paper_fp,
                                              offsite_dir = .get_config_path("LOCAL_SCEPTRE2_DATA_DIR"))
    modalities <- names(mm_odm@modalities)
    remaining_modalities <- modalities[!(modalities %in% c("grna_assignment", "grna_expression"))]
    grna_assign_modality <- get_modality(mm_odm, "grna_assignment")

    # load the gRNA assignments
    cellwise_grna_assignments <- grna_assign_modality |>
      get_cell_covariates() |>
      dplyr::pull(assigned_grna) 
    
    cellwise_grna_assignments_2 <- lowmoi:::get_grna_assignments_via_max_op(grna_assign_modality)      
    
    # NOT A SUBSET OF ROW NAMES OF FEATURE COVARIATES OF GRNA ASSIGN MODALITY
    unique_grnas <- unique(cellwise_grna_assignments)
    grna_tbl <- lapply(unique_grnas, function(grna) {
      which(cellwise_grna_assignments == grna)
    }) |> stats::setNames(nm = unique_grnas)
    n_cells <- sapply(grna_tbl, length)
    grna_target_df <- grna_assign_modality |> get_feature_covariates()
    grna_target_df <- data.frame(grna_id = factor(row.names(grna_target_df)),
                                 target = factor(grna_target_df$target))

    # loop through modalities
    x <- lapply(X = remaining_modalities, FUN = function(modality) {
      print(paste0("Working on modality ", modality))
      curr_modality <- get_modality(mm_odm, modality)
      feature_ids <- get_feature_ids(curr_modality)
      # loop through feature IDs; load the relevant feature
      lapply(X = seq(1, length(feature_ids)), FUN = function(i) {
      # lapply(X = seq(1, 4), FUN = function(i) {
        feature_id <- feature_ids[i]
        if (i %% 200 == 0) print(paste0("Working on ", feature_id, "."))
        gene_exp <- curr_modality[[feature_id,]] |> as.numeric()
        # loop through grnas, finding the number of nonzero cells for each grna
        n_nonzero_cells <- sapply(grna_tbl, FUN = function(curr_grna_idx) {
          (gene_exp[curr_grna_idx] >=1) |> sum()
        })
        df <- data.frame(feature_id = feature_id,
                         grna_id = names(grna_tbl),
                         n_nonzero_cells = n_nonzero_cells,
                         n_cells = n_cells)
        rownames(df) <- NULL
        return(df)
      }) |> data.table::rbindlist() |>
        dplyr::mutate(feature_id = factor(feature_id),
                      grna_id = factor(grna_id),
                      modality = factor(modality))
      }) |> data.table::rbindlist() |>
      dplyr::mutate(dataset = factor(dataset))
    y <- dplyr::left_join(x, grna_target_df, by = "grna_id") |>
      dplyr::mutate(target = factor(target))
  }) |> data.table::rbindlist() |>
    dplyr::mutate(paper = factor(paper))
}) |> data.table::rbindlist()

to_save_fp <- paste0(sceptre2_sample_sizes_dir, "n_nonzero_cells_per_grna.rds")
dataset_concat <- paste0(df$paper, "/", df$dataset, "/", df$modality) |> factor()
df <- df |> dplyr::mutate(dataset_concat = dataset_concat)
saveRDS(object = df, file = to_save_fp)
