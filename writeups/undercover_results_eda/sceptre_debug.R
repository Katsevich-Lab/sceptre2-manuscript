library(tidyverse)
load_all("~/research_code/lowmoi/")
load_all("~/research_code/sceptre2/")
sceptre_debug <- readRDS("/Users/timbarry/research_offsite/projects/sceptre2/results/undercover_grna_analysis/sceptre_debug.rds")

# fit levels
# excellent: < 1e-3
# good: > 1e-3, < 1e-2
# adequate: > 1e-2 < 5e-2
# poor: > 5e-2

# plot a few histograms from Frangieh co culture (gene)
pair <- sceptre_debug |>
  dplyr::filter(dataset == "frangieh/co_culture/gene",
                ks_stat > 1e-2, ks_stat < 5e-2) |>
  dplyr::sample_n(1) |>
  dplyr::select(response_id, undercover_grna)
response_id <- as.character(pair$response_id)
undercover_grna <- as.character(pair$undercover_grna)
sceptre_args <- lowmoi::get_sceptre_function_args_for_pair(response_id = response_id,
                                                           undercover_grna = undercover_grna,
                                                           dataset_name = "frangieh/co_culture/gene",
                                                           output_amount = 3,
                                                           B = 1000)
out <- do.call(what = sceptre2::run_sceptre_low_moi, args = sceptre_args)
sceptre2:::plot_fitted_density_result_row(out)



# first, plot the ks stats
ggplot(sceptre_debug, mapping = aes(x = ks_stat)) +
  geom_histogram() +
  facet_wrap(.~dataset) 
mean(sceptre_debug$ks_stat < 0.05) # 75\% of pairs have a KS stat of below 0.06, which seems good.

# focus on those pairs with good ks stats
sceptre_debug_good <- sceptre_debug |> filter(ks_stat < 0.05)
ggplot(data = sceptre_debug_good,
       mapping = aes(x = n_treatment_cells_with_expression, y = ks_stat)) +
  facet_wrap(.~dataset) + geom_point(alpha = 0.8, cex = 0.8) + scale_x_continuous(trans = "log10")
ggplot(data = sceptre_debug,
       mapping = aes(x = n_control_cells_with_expression, y = ks_stat)) +
  facet_wrap(.~dataset) + geom_point(alpha = 0.8, cex = 0.8) + scale_x_continuous(trans = "log10")

