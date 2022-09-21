# for each dataset, comptue the following:
# i) for each gRNA, the number of cells receiving that gRNA, as well as the number of cells with nonzero gene expression
# i) furthermore, the type of each gRNA
# Thus, we seek to create a data frame with the following columns:
# i) dataset, 
# ii) gRNA,
# iii) gene
# iii) n cells with gRNA
# iv) n with gRNA with nonzero gene expression
# v) gRNA type

library(ondisc)
sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")

papers <- c("frangieh", "liscovitch",  "papalexi", "schraivogel", "simulated")
for (paper in papers) {
  paper_dir <- paste0(sceptre2_data_dir, paper, "/")
  datasets <- list.files(paper_dir)
  for (dataset in datasets) {
    print(paste0("paper: ", paper, " dataset: ", dataset))
        
  }
}
