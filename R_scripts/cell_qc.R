library(ondisc)
sceptre2_offsite_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

# This script performs cell-wise QC.
# We perform cell-wise QC on datasets from two papers: Frangieh and Liscovitch.
# We follow as closely as possible the cell QC implemented by the original authors.
# We do not perform cell QC on the Schraivogel or Papalexi datasets, as these datasets already come with cell QC.

# 0) We write a general read and save multimodal_odm functions
# save
save_multimodal_odm <- function(multimodal_odm, paper, dataset, metadata_file_name, sceptre2_data_dir = paste0(sceptre2_offsite_dir, "data/")) {
  dataset_dir <- paste0(sceptre2_data_dir, paper, "/", dataset)
  modality_list <- multimodal_odm@modalities |> names()
  for (modality in modality_list) {
    save_odm(odm = get_modality(multimodal_odm, modality),
             metadata_fp = paste0(dataset_dir, "/", modality, "/", metadata_file_name))
  }
}

# read
read_multimodal_odm <- function(paper, dataset, sceptre2_data_dir = paste0(sceptre2_offsite_dir, "data/")) {
  dataset_dir <- paste0(sceptre2_data_dir, paper, "/", dataset)
  modality_vect <- list.files(dataset_dir)
  odm_list <- list()
  for (modality in modality_vect) {
    odm_dir <- paste0(dataset_dir, "/", modality, "/")
    curr_odm <- read_odm(odm_fp = paste0(odm_dir, "matrix.odm"),
                         metadata_fp = paste0(odm_dir, "metadata_orig.rds"))
    odm_list <- c(odm_list, curr_odm)
  }
  names(odm_list) <- modality_vect
  ret <- multimodal_ondisc_matrix(odm_list)
  return(ret)
}

# 1) Liscovitch
# First, load the experiment_small and experiment_big multimodal odms
exp_small_multimodal <- read_multimodal_odm(paper = "liscovitch", dataset = "experiment_small")
exp_big_multimodal <- read_multimodal_odm(paper = "liscovitch", dataset = "experiment_big")

# Second, define the QC function
run_lisc_cell_qc <- function(multimodal_odm, n_gRNA_reads = 100, n_atac_frags = 500, p_reads_one_gRNA_seq) {
  gRNA_mod <- get_modality(multimodal_odm, "grna")
  assignable <- gRNA_mod[[1:nrow(gRNA_mod),]] |> apply(MARGIN = 2, FUN = function(col) {
    max(col)/sum(col) > p_reads_one_gRNA_seq
  })
  multimodal_odm <- multimodal_odm[,assignable]
  high_quality <-  multimodal_odm |>
    get_cell_covariates() |> 
    dplyr::filter(chromatin_n_fragments >= n_atac_frags) |>
    dplyr::filter(grna_n_umis >= 100) |> row.names()
  multimodal_odm[,high_quality]
}

# Third, apply the QC, first to the small experiment and then to the big experiment
experiment_small_multimodal_cell_qc <- run_lisc_cell_qc(exp_small_multimodal, p_reads_one_gRNA_seq = 0.99)
experiment_big_multimodal_cell_qc <- run_lisc_cell_qc(exp_big_multimodal, p_reads_one_gRNA_seq = 0.90)

# Finally, save the multimodal odms
save_multimodal_odm(experiment_small_multimodal_cell_qc, "liscovitch", "experiment_small", "metadata_cell_qc.rds")
save_multimodal_odm(experiment_big_multimodal_cell_qc, "liscovitch", "experiment_big", "metadata_cell_qc.rds")

# 2) Frangieh
# First, load the co_culture, control, and ifn_gamma multimodal_odms
multimodal_odm_list <- list(co_culture = read_multimodal_odm(paper = "frangieh", dataset = "co_culture"),
                            control = read_multimodal_odm(paper = "frangieh", dataset = "control"),
                            ifn_gamma = read_multimodal_odm(paper = "frangieh", dataset = "ifn_gamma"))

# Second, write the QC function; this function simply restricts attention to cells that have received exactly 1 gRNA.
run_frangieh_cell_qc <- function(multimodal_odm) {
  cellwise_grna_count <- multimodal_odm |>
    get_cell_covariates() |>
    dplyr::pull(grna_n_nonzero)
  ret <- multimodal_odm[,cellwise_grna_count == 1]
  return(ret)
}

# Third, apply the QC function to the data
multimodal_odm_list_qc <- lapply(X = multimodal_odm_list, FUN = run_frangieh_cell_qc)

# Finally, save the multimodal odms
exp_names <- names(multimodal_odm_list_qc)
for (exp_name in exp_names) {
  curr_multimodal_odm <- multimodal_odm_list_qc[[exp_name]]
  save_multimodal_odm(multimodal_odm = curr_multimodal_odm,
                      paper = "frangieh",
                      dataset = exp_name,
                      metadata_file_name = "metadata_cell_qc.rds")
}
