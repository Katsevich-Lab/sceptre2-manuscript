library(tidyverse)
library(ondisc)

# path to SCEPTRE2 files
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
# create the multimodal odm
gasp_data_dir= paste0(sceptre2_dir, "data/gasperini/at_scale/")
# i) multimodal metadata file
multimodal_metadata_fp= paste0(gasp_data_dir, "multimodal_metadata.rds")
# ii) gene ODM
gene_odm_fp = paste0(gasp_data_dir, "gene/matrix.odm")
# iii) grna ODM
grna_odm_fp = paste0(gasp_data_dir, "grna_expression/matrix.odm")
# read the Gasperini data
mm_odm <- read_multimodal_odm(c(gene_odm_fp, grna_odm_fp), multimodal_metadata_fp)
gene_odm <- mm_odm |> get_modality("gene")
grna_odm <- mm_odm |> get_modality("grna_assignment")

cis_pairs <- read_tsv(paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/raw/GSE120861_gene_gRNAgroup_pair_table.at_scale.txt"))
cis_pairs |> 
  filter(general_group %in% c("TSS", "positive_ctrl", "NTC"),
         ENSG.targetgene %in% get_feature_ids(gene_odm),
         gRNAgroup %in% (get_feature_covariates(grna_odm) |> pull(grna_group))) |>
  rename(grna_group = gRNAgroup, gene_id = ENSG.targetgene) |>
  select(grna_group, gene_id, general_group) |>
  saveRDS(paste0(sceptre2_dir, "data/gasperini/at_scale/benchmarking_pairs.rds"))