source(paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R"))
library(tidyverse)
library(katlabutils)
library(cowplot)

full_base_dir <- "~/sceptre_outputs_lowmoi_full_stat"
calib_base_dir <- "~/sceptre_outputs_lowmoi_resid_stat"

# define several functions
join_result_dfs <- function(full_stat_df, resid_stat_df) {
  left_join(full_stat_df |> select(response_id, grna_target, p_value),
            resid_stat_df |> select(response_id, grna_target, p_value),
            by = c("response_id", "grna_target"),
            suffix = c("_full_stat", "_resid_stat"))
}

obtain_combined_df <- function(analysis_type, full_base_dir, calib_base_dir) {
  full_stat_results <- readRDS(paste0(full_base_dir, "/results_", analysis_type, ".rds"))
  resid_stat_results <- readRDS(paste0(calib_base_dir, "/results_", analysis_type, ".rds"))
  join_result_dfs(full_stat_results, resid_stat_results)
}

clip_result_df <- function(combined_discovery_analysis_df) {
  combined_discovery_analysis_df |>
    na.omit() |>
    mutate(p_value_full_stat = ifelse(p_value_full_stat < 1e-40, 1e-40, p_value_full_stat),
           p_value_resid_stat = ifelse(p_value_resid_stat < 1e-40, 1e-40, p_value_resid_stat))
}

# obtain the combined dfs
combined_calib_check_df <- obtain_combined_df("run_calibration_check", full_base_dir, calib_base_dir)
combined_power_check_df <- obtain_combined_df("run_power_check", full_base_dir, calib_base_dir)
combined_discovery_analysis_df <- obtain_combined_df("run_discovery_analysis", full_base_dir, calib_base_dir)

p1 <- ggplot(clip_result_df(combined_calib_check_df),
             aes(y = p_value_full_stat, x = p_value_resid_stat)) +
  geom_point(size = 0.9, col = "indianred2") +
  theme_bw() +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  geom_abline(slope = 1, intercept = 0) +
  ylab("p (score statistic)") +
  xlab("p (resid statistic)") +
  my_theme +
  ggtitle("Calibration check")

p2 <- ggplot(clip_result_df(combined_power_check_df),
             aes(y = p_value_full_stat, x = p_value_resid_stat)) +
  geom_point(col = "mediumseagreen") +
  theme_bw() +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  geom_abline(slope = 1, intercept = 0) +
  ylab("p (score statistic)") +
  xlab("p (resid statistic)") +
  my_theme +
  ggtitle("Positive control analysis")

p3 <- ggplot(data = clip_result_df(combined_discovery_analysis_df),
             aes(x = p_value_resid_stat, y = p_value_full_stat)) +
  geom_point(size = 0.9, col = "dodgerblue3") +
  theme_bw() +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("p (resid statistic)") +
  ylab("p (score statistic)") +
  my_theme +
  ggtitle("Discovery analysis")

p4 <- ggplot() + theme_minimal()
p_final <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = "auto")
