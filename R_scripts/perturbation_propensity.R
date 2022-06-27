data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- list.files(data_dir)

# loop over papers
results <- list()
for(paper in papers){
  results$paper <- list()
  paper_dir <- paste0(data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  # loop over datasets within a paper
  for (dataset in datasets) {
    results$paper$dataset <- list()

    # get gRNA ODM
    grna_fp <- paste0(paper, "/", dataset, "/grna_assignment")
    grna_odm <- lowmoi::load_dataset_modality(data_fp = gRNA_fp)

    # get gene ODM
    gene_fp <- paste0(paper, "/", dataset, "/gene")
    gene_odm <- lowmoi::load_dataset_modality(data_fp = gene_fp)

    # get list of NTC gRNAs
    ntcs <- gRNA_odm |>
      get_feature_covariates() |>
      filter(target == "non-targeting") |>
      rownames() |>
      unique()

    # get the gRNA assigned to each cell
    assigned_grnas <- sapply(1:ncol(gRNA_odm), function(col)(gRNA_odm[[,col]] |> as.logical() |> which()))
    grnas <- gRNA_odm |> get_feature_ids()
    assigned_grnas <- lapply(assigned_grnas, function(grna_idx)(if(length(grna_idx) == 1) grnas[grna_idx] else ""))


    # assigned_grnas <- gRNA_odm |>
    #   get_cell_covariates() |>
    #   pull(assigned_gRNA)


    # get a logical vector of whether a cell was assigned an NTC
    ntc_cells <- assigned_grnas |>
      mutate(ntc = assigned_gRNA %in% ntcs) |>
      pull(ntc)

    # restrict assigned gRNAs to cells with an NTC
    assigned_grnas <- assigned_grnas[ntc_cells]

    # extract cell covariates from gene ODM, and restrict attention to cells with NTC
    cell_covariates <- gene_odm[,ntc_cells] |> get_cell_covariates()

    # carry out the regression for each NTC
    for(ntc in ntcs){
      perturbation_indicator <- assigned_grnas == ntc
      glm_fit <- glm(perturbation_indicator ~ cell_covariates, family = "binomial")
    }
  }
}
