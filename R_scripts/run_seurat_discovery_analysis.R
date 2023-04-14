args <- commandArgs(trailingOnly = TRUE)
dataset_name <- args[1] # "schraivogel/enhancer_screen_chr11/gene", "schraivogel/enhancer_screen_chr8/gene", "papalexi/eccite_screen/gene", "frangieh/control/gene"
LOCAL_SCEPTRE2_DATA_DIR <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

library(ondisc)
library(lowmoi)

# load the response ODM and gRNA ODM
response_odm <- lowmoi::load_dataset_modality(dataset_name)
grna_dataset_name <- lowmoi::get_grna_dataset_name(dataset_name, "assignment")
grna_odm <- lowmoi::load_dataset_modality(grna_dataset_name)

# generate the pairs to analyze
grna_groups_to_keep <- grna_odm |>
  ondisc::get_feature_covariates() |>
  dplyr::group_by(target) |>
  dplyr::summarize(count = sum(n_nonzero)) |>
  dplyr::filter(count >= 3) |>
  dplyr::pull(target)
grna_groups_to_keep <- grna_groups_to_keep[grna_groups_to_keep != "non-targeting"]
response_grna_group_pairs <- expand.grid(response_id = ondisc::get_feature_ids(response_odm),
                                         grna_group = grna_groups_to_keep)

# run the method
res <- seurat_de(response_odm = response_odm,
                 grna_odm = grna_odm,
                 response_grna_group_pairs = response_grna_group_pairs)

# save the result
f_name <- paste0("seurat_", gsub(pattern = "/", x = dataset_name, fixed = TRUE, replacement = "_"), "_res.rds")
result_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "results/discovery_analyses/")
saveRDS(object = res, file = paste0(result_dir, f_name))
