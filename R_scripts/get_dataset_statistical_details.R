# load packages
library(tidyverse)
library(ondisc)

# load the shared figure script
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)

# set the result and data directories; load PC and undercover results
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
pc_res <- readRDS(paste0(result_dir, "positive_control_analysis/pc_results_processed.rds")) |>
  filter(n_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_control >= N_NONZERO_CONTROL_CUTOFF) 
undercover_res <- readRDS(paste0(result_dir, "undercover_grna_analysis/undercover_result_grp_1_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF)

# first, the number of negative control pairs
undercover_res |>
  filter(method == "seurat_de") |>
  group_by(dataset) |>
  summarize(count = n())

# next, the number of positive control pairs
pc_res |>
  filter(method == "seurat_de") |>
  group_by(dataset) |>
  summarize(count = n())

# loop over datasets, computing sample size information
papers <- c("frangieh",  "papalexi", "schraivogel", "simulated")
df <- lapply(papers, function(paper) {
  print(paste0("paper: ", paper))
  paper_dir <- paste0(data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  lapply(X = datasets, FUN = function(dataset) {
    print(paste0("paper: ", paper, " dataset: ", dataset))
    paper_fp <- paste0(paper, "/", dataset)
    mm_odm <- lowmoi::load_dataset_multimodal(paper_fp = paper_fp,
                                              offsite_dir = .get_config_path("LOCAL_SCEPTRE2_DATA_DIR"))
    modalities <- names(mm_odm@modalities)
    grna_modality <- mm_odm |> get_modality("grna_assignment")

    
    # get gRNA info
    grna_info <- grna_modality |> get_feature_covariates()
    n_nt_grnas <- grna_info |> filter(target_type == "non-targeting") |> nrow()
    n_targeting_grnas <- nrow(grna_info) - n_nt_grnas
    n_targeted_sites <- grna_info |> filter(target_type != "non-targeting") |> pull(target) |> unique() |> length()
    n_cells <- ncol(grna_modality)
    remaining_modalities <- modalities[!(modalities %in% c("grna_assignment", "grna_expression"))]
    
    # loop over remaining modalitites
    lapply(remaining_modalities, function(remaining_modality) {
      response_modality <- mm_odm |> get_modality(remaining_modality)
      # get response info
      n_responses <- nrow(response_modality)
      data.frame(n_responses = n_responses,
                 n_cells = n_cells,
                 n_targeting_grnas = n_targeting_grnas,
                 n_targeted_sites = n_targeted_sites,
                 n_nt_grnas = n_nt_grnas,
                 paper_fp = paper_fp,
                 modality = remaining_modality)
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist()
}) |> data.table::rbindlist()

df[c(1,3,5,7,8,9,10,13),]
