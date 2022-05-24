# This script performs a lightweight feature QC on all datasets under consideration.

# set directories
sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- list.files(sceptre2_data_dir)
# remove simulated paper
papers <- papers[papers != "simulated"]
# set params
N_CELLS_PER_GRNA_THRESH <- 10
FRAC_EXPRESSED_TRHESH <- 0.005

# load packages
library(ondisc)

# loop over papers
for (paper in papers) {
  # loop over datasets
  paper_dir <- paste0(sceptre2_data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  for (dataset in datasets) {
    # loop over modalities
    dataset_dir <- paste0(paper_dir, dataset, "/")
    modalities <- list.files(dataset_dir)
    for (modality in modalities) {
      print(paste0("paper: ", paper, ", dataset: ", dataset, ", modality: ", modality))
      modality_dir <- paste0(dataset_dir, modality, "/")
      metadata_fp <- paste0(modality_dir, "metadata_cell_qc.rds")
      if (!file.exists(metadata_fp)) metadata_fp <- paste0(modality_dir, "metadata_orig.rds")
      to_save_metadata_fp <- paste0(modality_dir, "metadata_qc.rds")
      odm_fp <- paste0(modality_dir, "matrix.odm")
      odm <- read_odm(odm_fp = odm_fp, metadata_fp = metadata_fp)
      odm_m <- lowmoi::load_whole_odm(odm)
      # if the modality is NOT gRNA...
      if (modality != "grna") {
        # check that the features have no underscores; if so, replace with dashes (both feature IDs and row names of feature covariate matrix)
        if (any(grepl(pattern = "_", x = get_feature_ids(odm), fixed = TRUE))) {
          odm@ondisc_matrix@feature_ids <- gsub(pattern = "_", replacement = "-",
                                                x = odm@ondisc_matrix@feature_ids, fixed = TRUE)
          row.names(odm@feature_covariates) <- gsub(pattern = "_", replacement = "-",
                                                    x = row.names(odm@feature_covariates), fixed = TRUE)
        }
        p_expressed <- Matrix::rowSums(odm_m >= 1)/ncol(odm_m)
        highly_exp_feats <- p_expressed >= FRAC_EXPRESSED_TRHESH
        # create a new metadata_fp if subset necessary
        odm_sub <- odm[highly_exp_feats,]
        save_odm(odm = odm_sub, metadata_fp = to_save_metadata_fp)
      } else {
        # the modality IS gRNA; filter out gRNAs expressed in fewer than 10 cells
        grna_mat <- lowmoi::load_whole_odm(odm)
        grna_assignments <- apply(X = grna_mat, MARGIN = 2, FUN = function(col) names(which.max(col)))
        grna_assignment_counts <- table(grna_assignments)
        good_grnas <- names(grna_assignment_counts[grna_assignment_counts >= N_CELLS_PER_GRNA_THRESH])
        odm_sub <- odm[good_grnas,]
        save_odm(odm = odm_sub, metadata_fp = to_save_metadata_fp)
      }
    }
  }
}
