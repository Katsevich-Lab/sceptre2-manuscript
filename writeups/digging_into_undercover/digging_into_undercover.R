library(ondisc)
library(lowmoi)
library(tidyverse)

res_fps <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                  "results/undercover_grna_analysis/undercover_result_grp_size_", 1:3,".rds")

# make histograms for Seurat and perm on all datasets and for all group sizes.
undercover_res_grp_1 <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                               "results/undercover_grna_analysis/undercover_result_grp_size_1.rds") |>
  readRDS() |>
  tibble::as_tibble() |>
  select(-clock_time, -max_ram) |>
  mutate(dataset_sub = sub(pattern = "_", replacement = "/", x = dataset, fixed = TRUE) |>
           stringi::stri_replace_last_fixed(str = _, pattern = "_", replacement = "/"))

undercover_res_sub <- undercover_res_grp_1 |> filter(method %in% c("seurat_de", "permutation_test"),
                                                     dataset %in% c("frangieh_co_culture_gene", "frangieh_ifn_gamma_gene", "frangieh_control_gene"))
undercover_res_sub_wide <- pivot_wider(data = undercover_res_sub,
                                       names_from = "method",
                                       values_from = p_value)

# make a simple histogram
ggplot(data = undercover_res_sub |> sample_frac(size = 0.05),
       mapping = aes(x = p_value)) +
  geom_histogram(bins = 20, col = "white", fill = "black") +
  facet_grid(method ~ dataset) +
  theme_bw()

# correlate the bulk p-values
blue_int <- 0.025
blue_slope <- 1.22
undercover_res_sub_wide_samp <- undercover_res_sub_wide |> sample_frac(size = 0.01)

ggplot(data = undercover_res_sub_wide_samp,
       mapping = aes(x = seurat_de, y = permutation_test)) +
  facet_wrap(dataset~.) +
  geom_point(size = 0.8, alpha = 0.3) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1,
              col = "darkred", lwd = 1.1) + 
  geom_abline(intercept = blue_int, slope = blue_slope,
              col = "blue", lwd = 1.1)

# filter the pairs that lie along the blue line
undercover_res_sub_wide_samp_near_line <- undercover_res_sub_wide_samp |>
  mutate(d_from_line = abs(blue_slope * seurat_de - permutation_test + blue_int)/sqrt(1 + blue_slope^2)) |>
  filter(d_from_line < 0.05)
ggplot(data = undercover_res_sub_wide_samp_near_line,
       mapping = aes(x = seurat_de, y = permutation_test)) + 
  facet_wrap(dataset~.) +
  geom_point(size = 0.8, alpha = 0.3) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1,
              col = "darkred", lwd = 1.1) # +
#  scale_x_continuous(trans = revlog_trans(base = 10)) +
#  scale_y_continuous(trans = revlog_trans(base = 10))
#  geom_abline(intercept = blue_int, slope = blue_slope,
#              col = "blue", lwd = 1.1) 

# pick a random point(s) to investigate
pair_to_check <- undercover_res_sub_wide_samp_near_line |>
   filter(seurat_de >= 0.45 & seurat_de <= 0.55) |>
   sample_n(1)


# load data
dataset_name <- pair_to_check$dataset_sub
undercover_ntc_name_in <- as.character(pair_to_check$undercover_grna)

if (FALSE) {
  dataset_name <- "frangieh/co_culture/gene"
  undercover_ntc_name <- "NO-SITE-812"
}

response_odm <- load_dataset_modality(dataset_name)
grna_dataset_name <- get_grna_dataset_name(dataset_name, "assignment")
grna_odm <- load_dataset_modality(grna_dataset_name)
undercover_ntc_name <- strsplit(x = undercover_ntc_name_in, split = ",", fixed = TRUE) |> unlist()
grna_feature_covariates <- grna_odm |> get_feature_covariates()
grna_feature_covariates[undercover_ntc_name, "target"] <- "undercover"
grna_odm_swapped <- grna_odm |> mutate_feature_covariates(target = grna_feature_covariates$target)
pairs_df <- pair_to_check |> summarize(grna_group = "undercover", response_id = pair_to_check$response_id)

# run permutation test
perm_res <- permutation_test(response_odm = response_odm,
                 grna_odm = grna_odm_swapped,
                 response_grna_group_pairs = pairs_df,
                 return_permuted_test_stats = TRUE,
                 test_stat = "mann_whit")
# run seurat
seurat_de(response_odm = response_odm,
          grna_odm = grna_odm_swapped,
          response_grna_group_pairs = pairs_df)


# plot resampled test stats
x <-perm_res |>
  select(resample_1:resample_1000) |>
  as.numeric()

