# load libraries and resolve conflicts
library(ondisc) 
library(sceptre)
library(readr)
library(dplyr)
library(conflicted)
conflicted::conflicts_prefer(dplyr::filter)

# set up directories
LOCAL_SCEPTRE2_DATA_DIR <-.get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
LOCAL_SCHRAIVOGEL_DATA_DIR <-.get_config_path("LOCAL_SCHRAIVOGEL_2020_DATA_DIR")
schraivogel_chr8_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, 
                               "data/schraivogel/enhancer_screen_chr8/")

# gene info
gene_odm_fp <- paste0(schraivogel_chr8_dir, "gene/matrix.odm")
gene_metadata_fp <- paste0(schraivogel_chr8_dir, "gene/metadata_qc.rds")
gene_odm <- read_odm(odm_fp = gene_odm_fp, metadata_fp = gene_metadata_fp)
gene_covariate_matrix <- gene_odm |> get_cell_covariates() 
gene_expression_matrix <- gene_odm[[seq(1, nrow(gene_odm)),]]
rownames(gene_expression_matrix) <- get_feature_ids(gene_odm)

# grna info
grna_odm_fp <- paste0(schraivogel_chr8_dir, "grna_expression/matrix.odm")
grna_metadata_fp <- paste0(schraivogel_chr8_dir, "grna_expression/metadata_qc.rds")
grna_odm <- read_odm(odm_fp = grna_odm_fp, metadata_fp = grna_metadata_fp)
grna_matrix <- grna_odm[[seq(1, nrow(grna_odm)),]]
grna_groups <- data.frame(grna_id = rownames(grna_odm@feature_covariates),
                          grna_group = grna_odm@feature_covariates$target)

# get Schraivogel results for chromosome 8 pairs
replace_periods <- function(str){
  stringr::str_replace(stringr::str_replace(str, "[.]", ":"), "[.]", "-")  
}
add_suffixes <- function(str){
  ifelse(str %in% c("CCNE2", "CPQ", "DSCC1", "FAM83A", "LRRCC1", 
                    "OXR1", "PHF20L1", "RIPK2", "STK3", "UBR5"),
         paste0(str, "-TSS"),
         ifelse(str %in% c("GATA1", "HS2", "MYC", "ZFPM2"),
                paste0(str, "-enh"),
                str))
}
schraivogel_results_fp <- paste0(LOCAL_SCHRAIVOGEL_DATA_DIR, 
                                 "raw/ftp/diff_expr_screen_nGenesCovar.csv")
schraivogel_results <- read_csv(schraivogel_results_fp) |>
  filter(sample == "8iScreen1",
         gene %in% get_feature_ids(gene_odm)) |> # genes that passed our QC
  mutate(grna_group = add_suffixes(replace_periods(perturbation))) |>
  rename(response_id = gene) |>
  select(response_id, grna_group, pvalue, logFC)

# set arguments for SCEPTRE
response_matrix <- gene_expression_matrix
grna_matrix <- grna_matrix
rownames(grna_matrix) <- ondisc::get_feature_ids(grna_odm)
covariate_data_frame <- gene_covariate_matrix
grna_group_data_frame <- grna_groups
formula_object <- ~log(n_umis) + log(n_nonzero) # + batch
calibration_check <- FALSE
response_grna_group_pairs <- schraivogel_results |> select(response_id, grna_group)

# run SCEPTRE
result_sceptre <- run_sceptre_lowmoi(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  covariate_data_frame = covariate_data_frame,
  grna_group_data_frame = grna_group_data_frame,
  formula_object = formula_object,
  response_grna_group_pairs = response_grna_group_pairs,
  calibration_check = calibration_check, 
  return_debugging_metrics = TRUE
)

# save the results
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
output_dir <- paste0(sceptre2_dir, "results/schraivogel_analysis/")
if (!dir.exists(output_dir)) dir.create(output_dir)
saveRDS(result_sceptre, paste0(output_dir, "sceptre_schraivogel_chr_8_results.rds"))
saveRDS(schraivogel_results, paste0(output_dir, "schraivogel_schraivogel_chr_8_results.rds"))