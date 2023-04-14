args <- commandArgs(trailingOnly = TRUE)
dataset_name <- args[1] # "schraivogel/enhancer_screen_chr11/gene", "schraivogel/enhancer_screen_chr8/gene", "papalexi/eccite_screen/gene", "frangieh/control/gene"
LOCAL_SCEPTRE2_DATA_DIR <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

library(ondisc)
library(lowmoi)

response_odm <- lowmoi::load_dataset_modality(dataset_name)
grna_dataset_name <- lowmoi::get_grna_dataset_name(dataset_name, "assignment")
grna_odm <- lowmoi::load_dataset_modality(grna_dataset_name)
response_grna_group_pairs <- generate_all_pairs(response_odm, grna_odm)
res <- seurat_de(response_odm = response_odm,
                 grna_odm = grna_odm,
                 response_grna_group_pairs = response_grna_group_pairs)

# save result
f_name <- paste0("seurat_", gsub(pattern = "/", x = dataset_name, fixed = TRUE, replacement = "_"), "_res.rds")
result_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "results/discovery_analyses/")
saveRDS(object = res, file = paste0(result_dir, f_name))
