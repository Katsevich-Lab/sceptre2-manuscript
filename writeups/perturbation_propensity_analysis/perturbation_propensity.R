data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- list.files(data_dir)

# set up results tibble
results <- tibble::tibble(
  paper = character(),
  dataset = character(),
  ntc = character(),
  covariate = character(),
  test_type = character(),
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

    # get data frame containing the perturbation indicators and covariates
    df <- lowmoi::get_data_for_pert_prop(paper, dataset)

    # get list of NTCs
    ntcs <- df |> dplyr::pull(assigned_grna) |> unique()

    # carry out the regression for each NTC
    for (ntc in ntcs) {
      # create perturbation indicator based on given NTC
      df_ntc <- df |>
        dplyr::mutate(perturbation_indicator = as.numeric(assigned_grna == ntc))

      # remove unnecessary columns from cell covariates
      df_ntc <- df_ntc |> dplyr::select(-assigned_grna)
      if (paper == "frangieh") {
        df_ntc <- df_ntc |> dplyr::select(-condition)
      }

      ### joint analysis ###

      # fit joint GLM
      glm_fit <- glm(perturbation_indicator ~ ., family = "binomial", data = df_ntc)

      # fit null GLM
      glm_fit_null <- glm(perturbation_indicator ~ 1, family = "binomial", data = df_ntc)

      # extract joint fit information
      fit_info <- coef(summary(glm_fit)) |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "covariate") |>
        dplyr::filter(covariate != "(Intercept)") |>
        dplyr::rename(estimate = Estimate,
               std_error = `Std. Error`,
               zvalue = `z value`,
               pvalue = `Pr(>|z|)`) |>
        dplyr::mutate(paper = paper,
               dataset = dataset,
               ntc = ntc,
               test_type = "joint") |>
        dplyr::select(paper, dataset, ntc, covariate, test_type,
               estimate, std_error, zvalue, pvalue)

      # add to results data frame
      results <- results |> dplyr::bind_rows(fit_info)

      # extract p-value for chi-squared test
      chisq_pval <- anova(glm_fit, glm_fit_null, test = "Chisq")[2, "Pr(>Chi)"]
      results <- results |> dplyr::add_row(paper = paper, dataset = dataset, ntc = ntc,
                                    covariate = "all", test_type = "joint", estimate = NA,
                                    std_error = NA, zvalue = NA, pvalue = chisq_pval)

      ### marginal analysis ###
      covariates <- df_ntc |> colnames() |> setdiff("perturbation_indicator")
      for(covariate in covariates){
        # fit marginal GLM
        glm_fit <- glm(sprintf("perturbation_indicator ~ %s", covariate),
                       family = "binomial", data = df_ntc)

        # extract marginal fit information
        fit_info <- coef(summary(glm_fit)) |>
          as.data.frame() |>
          tibble::rownames_to_column(var = "covariate") |>
          dplyr::filter(covariate != "(Intercept)") |>
          dplyr::rename(estimate = Estimate,
                 std_error = `Std. Error`,
                 zvalue = `z value`,
                 pvalue = `Pr(>|z|)`) |>
          dplyr::mutate(paper = paper,
                 dataset = dataset,
                 ntc = ntc,
                 test_type = "marginal") |>
          dplyr::select(paper, dataset, ntc, covariate, test_type,
                 estimate, std_error, zvalue, pvalue)

        # add to results data frame
        results <- results |> dplyr::bind_rows(fit_info)
      }
    }
  }
}

# save results to disk
results_dir <- paste0(
  .get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
  "results/perturbation_propensity_analysis"
)
dir.create(results_dir)
saveRDS(results, paste0(results_dir, "/results.rds"))
