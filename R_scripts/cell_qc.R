library(ondisc)
sceptre2_offsite_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

# This script performs cell-wise QC.
# We perform cell-wise QC on datasets from two papers: Frangieh and Liscovitch.
# We follow as closely as possible the cell QC implemented by the original authors.
# We do not perform cell QC on the Schraivogel or Papalexi datasets, as these datasets already come with cell QC.

# 0) We write a general save multimodal_odm function.
save_multimodal_odm <- function(multimodal_odm, paper, dataset, sceptre2_data_dir = paste0(sceptre2_offsite_dir, "data/")) {
  dataset_dir <- paste0(sceptre2_data_dir, paper, "/", dataset)
  modality_list <- multimodal_odm@modalities |> names()
  for (modality in modality_list) {
    save_odm(odm = get_modality(multimodal_odm, modality),
             metadata_fp = paste0(dataset_dir, "/", modality, "/metadata_qc.rds"))
  }
}

# 1) Liscovitch
# Load the chromatin and grna modalities for experiment_small and experiment_big
liscovitch_dir <- .get_config_path("LOCAL_LISCOVITCH_2021_DATA_DIR")

experiment_small_chromatin <- read_odm(odm_fp = paste0(liscovitch_dir, "processed/experiment_small/chromatin/chip_counts.odm"),
                                       metadata_fp = paste0(liscovitch_dir, "processed/experiment_small/chromatin/metadata.rds"))
experiment_small_grna <- read_odm(odm_fp = paste0(liscovitch_dir, "processed/experiment_small/gRNA/gRNA_counts.odm"),
                                  metadata_fp = paste0(liscovitch_dir, "processed/experiment_small/gRNA/metadata.rds"))
experiment_small_multimodal <- multimodal_ondisc_matrix(list(chromatin = experiment_small_chromatin, grna = experiment_small_grna))
experiment_big_chromatin <- read_odm(odm_fp = paste0(liscovitch_dir, "processed/experiment_big/chromatin/chip_counts.odm"),
                                       metadata_fp = paste0(liscovitch_dir, "processed/experiment_big/chromatin/metadata.rds"))
experiment_big_grna <- read_odm(odm_fp = paste0(liscovitch_dir, "processed/experiment_big/gRNA/gRNA_counts.odm"),
                                  metadata_fp = paste0(liscovitch_dir, "processed/experiment_big/gRNA/metadata.rds"))
experiment_big_multimodal <- multimodal_ondisc_matrix(list(chromatin = experiment_big_chromatin, grna = experiment_big_grna))

# define the general QC function
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

# small experiment cell QC (QC metrics: (a) >= 100 gRNA reads, (b) >= 99% of reads assigned to one gRNA sequence, (c) >= 500 ATAC-seq fragments)
experiment_small_multimodal_cell_qc <- run_lisc_cell_qc(experiment_small_multimodal, p_reads_one_gRNA_seq = 0.99)
# big experiment cell QC (QC metrics: >= 100 gRNA reads, (b) >= 90% of reads assigned to one gRNA sequence, (c) >= 500 ATAC-seq fragments)
experiment_big_multimodal_cell_qc <- run_lisc_cell_qc(experiment_big_multimodal, p_reads_one_gRNA_seq = 0.90)
# save multimodal odms
save_multimodal_odm(experiment_small_multimodal_cell_qc, "liscovitch", "experiment_small")
save_multimodal_odm(experiment_big_multimodal_cell_qc, "liscovitch", "experiment_big")

