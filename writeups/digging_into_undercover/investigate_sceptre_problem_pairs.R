fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/writeups/digging_into_undercover/figs_aug_2022/")
load_all("~/research_code/sceptre2/")
problem_pairs <- readRDS(file = "~/Desktop/sceptre_problem_pairs_frangieh_control_gene.rds")
control_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/frangieh/control/")
undercover_grnas <- problem_pairs$undercover_grna |> unique()

get_sceptre_function_args_for_pair <- function(problem_pairs, undercover_grna) {
  problem_responses <- problem_pairs |>
    dplyr::filter(undercover_grna == !!undercover_grna) |>
    dplyr::pull(response_id) |>
    as.character()

  dataset_name <- "frangieh/control/gene"
  undercover_ntc_name_in <- undercover_grna
  grna_modality <- "assignment"

  response_odm <- lowmoi::load_dataset_modality(dataset_name)
  grna_dataset_name <- lowmoi::get_grna_dataset_name(dataset_name, grna_modality)
  grna_odm <- lowmoi::load_dataset_modality(grna_dataset_name)

  undercover_ntc_name <- strsplit(x = undercover_ntc_name_in, split = ",", fixed = TRUE) |> unlist()
  grna_feature_covariates <- grna_odm |> ondisc::get_feature_covariates()
  grna_feature_covariates[undercover_ntc_name, "target"] <- "undercover"
  if (!("non-targeting" %in% grna_feature_covariates$target)) {
    stop("After performing label swap, `non-targeting` is no longer string in the `target` column.")
  }
  grna_odm_swapped <- grna_odm |> ondisc::mutate_feature_covariates(target = grna_feature_covariates$target)
  response_grna_group_pairs <- data.frame(response_id = ondisc::get_feature_ids(response_odm),
                                          grna_group = "undercover") |>
    dplyr::filter(response_id %in% problem_responses)
  response_odm <- response_odm; grna_odm <- grna_odm_swapped; response_grna_group_pairs <- response_grna_group_pairs
  mm_odm <- lowmoi::process_multimodal_odm(ondisc::multimodal_ondisc_matrix(list(response = response_odm, grna = grna_odm)))
  form <- mm_odm@modalities$response@misc$sceptre_formula

  ret <- list(mm_odm = mm_odm,
              response_grna_group_pairs = response_grna_group_pairs,
              form = form,
              response_modality_name = "response",
              grna_modality_name = "grna",
              grna_group_column_name = "target",
              B = 1e5,
              side = "both",
              full_output = 2)
  return(ret)
}

################
# ANALYSIS START
################
load_all("~/research_code/sceptre2/")
curr_undercover_grna <- as.character(sample(undercover_grnas, 1))
funct_args <- get_sceptre_function_args_for_pair(problem_pairs = problem_pairs,
                                                 undercover_grna = curr_undercover_grna)
funct_args$response_grna_group_pairs <- funct_args$response_grna_group_pairs |> dplyr::sample_n(15)
curr_problem_pairs <- problem_pairs |>
  dplyr::filter(undercover_grna == curr_undercover_grna)
out <- do.call(what = run_sceptre_low_moi, args = funct_args)

### going inside the function
mm_odm <- funct_args[[1]]
response_grna_group_pairs <- funct_args[[2]]
form <- funct_args[[3]]
response_modality_name <- funct_args[[4]]
grna_modality_name <- funct_args[[5]]
grna_group_column_name <- funct_args[[6]]
B <- funct_args[[7]]
side <- funct_args[[8]]
full_output <- funct_args[[9]]
###


ggplot(data = out, mapping = aes(x = n_control_cells_with_expression, y = ks_stat)) +
  geom_point() +
  scale_x_log10()

ggplot(data = out, mapping = aes(x = n_treatment_cells_with_expression, y = ks_stat)) +
  geom_point() +
  scale_x_log10()
