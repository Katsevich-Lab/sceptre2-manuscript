# This script performs a lightweight feature QC on all datasets under consideration.

# set directories
sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- list.files(sceptre2_data_dir)
# remove simulated paper
papers <- papers[papers != "simulated"]

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
      # if the modality is NOT gRNA... (if it is, then skip; no feature QC on the gRNA matrix for now)
      if (modality != "grna") {
        odm_fp <- paste0(modality_dir, "matrix.odm")
        # read the odm
        odm <- read_odm(odm_fp = odm_fp, metadata_fp = metadata_fp)
        # check that the features have no underscores; if so, replace with dashes (both feature IDs and row names of feature covariate matrix)
        if (any(grepl(pattern = "_", x = get_feature_ids(odm), fixed = TRUE))) {
          odm@ondisc_matrix@feature_ids <- gsub(pattern = "_", replacement = "-",
                                                x = odm@ondisc_matrix@feature_ids, fixed = TRUE)
          row.names(odm@feature_covariates) <- gsub(pattern = "_", replacement = "-",
                                                    x = row.names(odm@feature_covariates), fixed = TRUE)
        }
        highly_exp_feats <- get_highly_expressed_features(odm, frac_expressed = 0.005)
        # create a new metadata_fp if subset necessary
        odm_sub <- odm[highly_exp_feats,]
        save_odm(odm = odm_sub, metadata_fp = to_save_metadata_fp)
      } else {
        # the modality IS gRNA; simply create a symbolic link (no qc)
        system(paste("ln -s", metadata_fp, to_save_metadata_fp))
      }
    }
  }
}
