library(ondisc)
library(lowmoi)
library(dplyr)
library(tibble)

data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- c("frangieh", "papalexi", "schraivogel")
methods <- c(schraivogel_method, liscovitch_method, mimosca,
             weissman_method, seurat_de, permutation_test, nb_regression)
names(methods) <- c("schraivogel", "liscovitch", "mimosca",
                    "weissman", "seurat", "permutation",
                    "nb_regression")

# set up results tibble
results <- tibble(
  paper = character(),
  dataset = character(),
  method = character(),
  response_id = character(),
  grna_group = character(),
  p_value = numeric()
)

for (paper in papers) {
  paper_dir <- paste0(data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  # loop over datasets within a paper
  for (dataset in datasets) {
    cat(sprintf("Working on the %s dataset from %s...\n", dataset, paper))

    grna_fp <- paste0(paper, "/", dataset, "/grna_assignment")
    grna_odm <- load_dataset_modality(data_fp = grna_fp)

    response_fp <- paste0(paper, "/", dataset, "/gene")
    response_odm <-load_dataset_modality(data_fp = response_fp)

    if(paper == "schraivogel"){
      grna_odm <- grna_odm |> mutate_feature_covariates(target = known_effect)
    }
    targeted_genes <- intersect(grna_odm |> get_feature_covariates() |> pull(target),
                                  response_odm |> get_feature_ids())

    response_grna_group_pairs <- grna_odm |>
      get_feature_covariates() |>
      rownames_to_column(var = "grna_group") |>
      filter(target %in% targeted_genes) |>
      select(grna_group, target) |>
      rename(response_id = target)

    feature_covariates <- grna_odm |>
      get_feature_covariates() |>
      rownames_to_column(var = "grna") |>
      mutate(target = ifelse(target == "non-targeting",
                             "non-targeting",
                             grna)) |>
      column_to_rownames(var = "grna")

    grna_odm <- grna_odm |>
      mutate_feature_covariates(feature_covariates)

    for(method_name in names(methods)){
      cat(sprintf("Applying %s...\n", method_name))
      result <- do.call(methods[[method_name]], list(response_odm, grna_odm, response_grna_group_pairs[1,]))
      result <- result |>
        mutate(method = method_name, dataset = dataset, paper = paper) |>
        select(paper, dataset, method, response_id, grna_group, p_value)
      results <- results |> bind_rows(result)
    }
  }
}

# save results to disk
results_dir <- paste0(
  .get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
  "results/positive_control_analysis"
)
dir.create(results_dir)
saveRDS(results, paste0(results_dir, "/results.rds"))
