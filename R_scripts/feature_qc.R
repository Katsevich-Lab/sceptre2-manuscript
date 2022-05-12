# This script performs a lightweight feature QC on all datasets under consideration.

# set directories
sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
papers <- list.files(sceptre2_data_dir)

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
        highly_exp_feats <- get_highly_expressed_features(odm, frac_expressed = 0.005)
        if (nrow(odm) == length(highly_exp_feats)) {
          # symbolic link to the current metadata_fp if no subset necessary
          system(paste("ln -s", metadata_fp, to_save_metadata_fp))
        } else {
          # create a new metadata_fp if subset necessary
          odm_sub <- odm[highly_exp_feats,]
          save_odm(odm = odm_sub, metadata_fp = to_save_metadata_fp)
        }
      } else {
        # the modality IS gRNA; simply create a symbolic link (no qc)
        system(paste("ln -s", metadata_fp, to_save_metadata_fp))
      }
    }
  }
}
