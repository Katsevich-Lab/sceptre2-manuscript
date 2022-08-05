library(ondisc)
library(lowmoi)

response_id <- "FOXF2"
undercover_ntc_name_in <- "ONE-NON-GENE-SITE-328"
dataset_name <- "frangieh/co_culture/gene"

response_odm <- load_dataset_modality(dataset_name)
grna_dataset_name <- get_grna_dataset_name(dataset_name, "assignment")
grna_odm <- load_dataset_modality(grna_dataset_name)
undercover_ntc_name <- strsplit(x = undercover_ntc_name_in, split = ",", fixed = TRUE) |> unlist()
grna_feature_covariates <- grna_odm |> get_feature_covariates()
grna_feature_covariates[undercover_ntc_name, "target"] <- "undercover"
grna_odm_swapped <- grna_odm |> mutate_feature_covariates(target = grna_feature_covariates$target)
pairs_df <- data.frame(grna_group = "undercover", response_id = response_id)

perm_res_1 <- permutation_test(response_odm = response_odm,
                               grna_odm = grna_odm_swapped,
                               response_grna_group_pairs = pairs_df,
                               return_permuted_test_stats = TRUE,
                               test_stat = "log_fold_change")

# perm_res_2 <- permutation_test_pl


perm_res[,1:5]
perm_res |> dplyr::select(resample_1:resample_1000) |> as.numeric() |> hist()
