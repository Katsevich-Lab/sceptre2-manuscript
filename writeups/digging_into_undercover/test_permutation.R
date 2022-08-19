# load data, functions, and packages
library(ondisc)
library(lowmoi)
funct_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/writeups/digging_into_undercover/analyze_undercover_results_plot_functs.R")
source(funct_script)
undercover_res_fps <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/undercover_grna_analysis/undercover_result_grp_size_", seq(1, 3), ".rds")
undercover_res <- readRDS(undercover_res_fps[1]) |>
  dplyr::mutate(clock_time = NULL, max_ram = NULL)|>
  dplyr::mutate(dataset_slash = replace_dataset_underscore_with_slash(dataset))


undercover_res |>
  dplyr::filter(dataset_slash == "frangieh/co_culture/gene", method == "permutation_test") |>
  dplyr::arrange(p_value) |> head(20)

response_id <- "C1orf159"
undercover_ntc_name_in <- "NO-SITE-812"
dataset_name <- "frangieh/co_culture/gene"

pieces <- get_undercover_ingredients(response_id, undercover_ntc_name_in, dataset_name)
response_odm <- pieces$response_odm
grna_odm_swapped <- pieces$grna_odm_swapped
pairs_df <- pieces$pairs_df

perm_res_1 <- permutation_test(response_odm = response_odm,
                               grna_odm = grna_odm_swapped,
                               response_grna_group_pairs = pairs_df,
                               return_permuted_test_stats = TRUE,
                               test_stat = "log_fold_change")
if (perm_res_1$z_value == -Inf) perm_res_1$z_value <- -3

perm_res_2 <- permutation_test_plus(response_odm = response_odm,
                                    grna_odm = grna_odm_swapped,
                                    response_grna_group_pairs = pairs_df,
                                    return_permuted_test_stats = TRUE,
                                    n_rep = 10000)

perm_res_1[,1:5]
perm_res_2[,1:5]

sceptre::plot_result(perm_res_1)
sceptre::plot_result(perm_res_2)
