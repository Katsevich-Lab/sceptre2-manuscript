data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- list.files(data_dir)

# set up results tibble
results <- tibble(
  paper = character(),
  dataset = character(),
  ntc = character(),
  covariate = character(),
  estimate = numeric(),
  std_error = numeric(),
  zvalue = numeric(),
  pvalue = numeric()
)

# loop over papers
for (paper in papers) {
  paper_dir <- paste0(data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  # loop over datasets within a paper
  for (dataset in datasets) {
    cat(sprintf("Working on the %s dataset from %s...\n", dataset, paper))

    # get gRNA ODM
    grna_fp <- paste0(paper, "/", dataset, "/grna_assignment")
    grna_odm <- lowmoi::load_dataset_modality(data_fp = grna_fp)

    # get response ODM
    if (paper == "liscovitch") {
      response_fp <- paste0(paper, "/", dataset, "/chromatin")
    } else {
      response_fp <- paste0(paper, "/", dataset, "/gene")
    }
    response_odm <- lowmoi::load_dataset_modality(data_fp = response_fp)

    # get list of NTC gRNAs
    ntcs <- grna_odm |>
      get_feature_covariates() |>
      filter(target %in% c("non-targeting", "NT")) |>
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
    assigned_grnas <- assigned_grnas_df[ntc_cells, , drop = F] |> pull(assigned_gRNA)

    # extract cell covariates from gene ODM, and restrict attention to cells with NTC
    cell_covariates <- response_odm[, ntc_cells] |> get_cell_covariates()

    # append assigned_grnas column to cell_covariates
    df <- cell_covariates |>
      mutate(assigned_grna = assigned_grnas)

    # carry out the regression for each NTC
    for (ntc in ntcs) {
      # create perturbation indicator based on given NTC
      df_ntc <- df |>
        mutate(perturbation_indicator = as.numeric(assigned_grna == ntc))

      # remove unnecessary columns from cell covariates
      df_ntc <- df_ntc |> select(-assigned_grna)
      if (paper == "frangieh") {
        df_ntc <- df_ntc |> select(-condition)
      }
      if (paper == "papalexi") {
        df_ntc <- df_ntc |> select(-guide_ID)
      }

      # fit GLM
      glm_fit <- glm(perturbation_indicator ~ ., family = "binomial", data = df_ntc)

      # extract fit information
      fit_info <- coef(summary(glm_fit)) |>
        as.data.frame() |>
        rownames_to_column(var = "covariate") |>
        filter(covariate != "(Intercept)") |>
        rename(estimate = Estimate, std_error = `Std. Error`, zvalue = `z value`, pvalue = `Pr(>|z|)`) |>
        mutate(paper = paper, dataset = dataset, ntc = ntc) |>
        select(paper, dataset, ntc, covariate, estimate, std_error, zvalue, pvalue)

      # add to results data frame
      results <- results |> bind_rows(fit_info)
    }
  }
}

# save results to disk
results_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                      "results/perturbation_propensity_analysis")
dir.create(results_dir)
saveRDS(results, paste0(results_dir, "/results.rds"))
