# compute table in figure 5b
discovery_analysis_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/discovery_analyses/")
papalexi_res <- readRDS(paste0(discovery_analysis_dir, "papalexi_gene_calibration_res.rds"))
frangieh_res <- readRDS(paste0(discovery_analysis_dir, "frangieh_control_calibration_res.rds"))

n_papa_pairs <- nrow(papalexi_res)
papa_discovery_time <- 6382.82
papa_calibration_time <- 1942.02

n_frangieh_pairs <- nrow(frangieh_res)
frangieh_discovery_time <- 31020.98
frangieh_calibration_time <- 2721.97

# papalexi
n_papa_pairs/papa_discovery_time
n_papa_pairs/papa_calibration_time

# frangieh
n_frangieh_pairs/frangieh_discovery_time
n_frangieh_pairs/frangieh_calibration_time
