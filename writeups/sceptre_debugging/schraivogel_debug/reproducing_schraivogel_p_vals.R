# load packages
library(ondisc)
library(lowmoi)

###########
# LOAD ODMS
###########
LOCAL_SCEPTRE2_DATA_DIR <-.get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
schraivogel_chr8_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR,
                               "data/schraivogel/enhancer_screen_chr8/")
# response info
response_odm_fp <- paste0(schraivogel_chr8_dir, "gene/matrix.odm")
response_metadata_fp <- paste0(schraivogel_chr8_dir, "gene/metadata_qc.rds")
response_odm <- read_odm(odm_fp = response_odm_fp, metadata_fp = response_metadata_fp)

# grna info
grna_odm_fp <- paste0(schraivogel_chr8_dir, "grna_expression/matrix.odm")
grna_metadata_fp <- paste0(schraivogel_chr8_dir, "grna_expression/metadata_qc.rds")
grna_odm <- read_odm(odm_fp = grna_odm_fp, metadata_fp = grna_metadata_fp)

###################################
# LOAD ORIGINAL SSCHRAIVOGEL RESULT
###################################
orig_schraivogel_res <- readRDS(paste0(schraivogel_chr8_result_dir, "schraivogel_schraivogel_chr_8_results.rds"))
pairs_to_analyze <- data.frame(response_id = "PTDSS1",
                               grna_group = "chr8:97842054-97842529")

schraivogel_res_tim <- schraivogel_method(response_odm = response_odm,
                                          grna_odm = grna_odm,
                                          response_grna_group_pairs = pairs_to_analyze)
schraivogel_res_tim

orig_schraivogel_res |> dplyr::filter(response_id == pairs_to_analyze$response_id,
                                      grna_group == pairs_to_analyze$grna_group)
