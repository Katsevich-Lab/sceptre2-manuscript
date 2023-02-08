library(ondisc)
sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
paper <- "gasperini"
dataset <- "at_scale"
FRAC_EXPRESSED_TRHESH <- 0.005

# load the dataset into a multimodal ODM
print(paste0("paper: ", paper, " dataset: ", dataset))
paper_dir <- paste0(sceptre2_data_dir, paper, "/")
multimodal_metadata_fp <- paste0(paper_dir, dataset, "/multimodal_metadata.rds")
if (file.exists(multimodal_metadata_fp)) file.remove(multimodal_metadata_fp)
mm_odm <- lowmoi::read_all_modalities(paper, dataset)

# 1. add global formula obj and MOI
mm_odm@global_misc[["moi"]] <- "high"
mm_odm@global_misc[["formula"]] <- ~ log(gene_n_umis) + log(grna_expression_n_umis) + p_mito + batch 

# 2. perform basic feature QC on gene modality
gene_modality <- mm_odm |> get_modality("gene")
feats_to_keep <- get_highly_expressed_features(gene_modality, FRAC_EXPRESSED_TRHESH)
mm_odm@modalities$gene <- gene_modality[feats_to_keep,]

# 3. perform cell QC -- remove cells with 0 gRNA UMIs
good_cells <- mm_odm |>
  ondisc::get_cell_covariates() |> 
  dplyr::filter(grna_expression_n_umis >= 1) |>
  row.names()
mm_odm <- mm_odm[,good_cells]

# 3. ensure that `grna_group` column contains "non-targeting"
updated_grna_odm <- mm_odm |>
  ondisc::get_modality("grna_expression") |>
  ondisc::mutate_feature_covariates(grna_group = ifelse(target_type == "non-targeting",
                                                        "non-targeting", grna_group))
mm_odm@modalities[["grna_expression"]] <- updated_grna_odm

# 4. create a multimodal ondisc matrix free of redundancy and write
mm_odm_sub_proc <- lowmoi::process_multimodal_odm(mm_odm)
save_multimodal_odm(multimodal_odm = mm_odm_sub_proc,
                    multimodal_metadata_fp = multimodal_metadata_fp)
