data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- list.files(data_dir)

# loop over papers
results <- tibble(paper = character(), dataset = character(), ntc = character(),
                  pvalue = numeric(), zvalue = numeric(), effect_size = numeric())
for(paper in papers){
  paper_dir <- paste0(data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  # loop over datasets within a paper
  for (dataset in datasets) {
    # get gRNA ODM
    grna_fp <- paste0(paper, "/", dataset, "/grna_assignment")
    grna_odm <- lowmoi::load_dataset_modality(data_fp = gRNA_fp)

    # get gene ODM
    gene_fp <- paste0(paper, "/", dataset, "/gene")
    gene_odm <- lowmoi::load_dataset_modality(data_fp = gene_fp)

    # get list of NTC gRNAs
    ntcs <- grna_odm |>
      get_feature_covariates() |>
      filter(target == "non-targeting") |>
      rownames() |>
      unique()

    # get the gRNA assigned to each cell
    assigned_grnas_df <- grna_odm |>
      get_cell_covariates() |>
      select(assigned_gRNA)

    # get a logical vector of whether a cell was assigned an NTC
    ntc_cells <- assigned_grnas_df |>
      mutate(ntc = assigned_gRNA %in% ntcs) |>
      pull(ntc)

    # restrict assigned gRNAs to cells with an NTC
    assigned_grnas <- assigned_grnas_df[ntc_cells,,drop=F] |> pull(assigned_gRNA)

    # extract cell covariates from gene ODM, and restrict attention to cells with NTC
    cell_covariates <- gene_odm[,ntc_cells] |> get_cell_covariates()

    # append assigned_grnas column to cell_covariates
    df <- cell_covariates |>
      mutate(assigned_grna = assigned_grnas)

    # carry out the regression for each NTC
    for(ntc in ntcs){
      df_ntc <- df |>
        mutate(perturbation_indicator = as.numeric(assigned_grna == ntc)) |>
        select(-assigned_grna, -condition)
      glm_fit <- glm(perturbation_indicator ~ ., family = "binomial", data = df_ntc)
      fit_info <- coef(summary(glm_fit))
      results |> add_row(paper = paper, dataset = dataset, ntc = ntc,
                         pvalue = numeric(), zvalue = numeric(), effect_size = numeric())
    }
  }
}
