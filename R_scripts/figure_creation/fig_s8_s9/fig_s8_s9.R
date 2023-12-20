####################################################################
# The purpose of figures s8 and s9 is to compare the score statistic
# against the residual test statistic on real and simulated data.
####################################################################
# load packages
library(tidyverse); conflicts_prefer(dplyr::filter)
library(katlabutils)
library(cowplot)

# set file paths, load results, source helpers
source(paste0(.get_config_path("LOCAL_CODE_DIR"),
              "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R"))
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

############################
# PLOT 1: SIMULATION RESULTS
############################
# load simulation result, do basic preprocessing
simulation_result <- paste0(sceptre2_dir, "results/extra_analyses/score_vs_resid_sim.rds") |> readRDS()
cols <- c("deepskyblue4", "firebrick2")
simulation_result <- simulation_result$result_df |>
  mutate(label = factor(x = null_true, levels = c(FALSE, TRUE),
                        labels = c("Alternative true", "Null true")))

fit_1 <- lm(p_lrt_trans ~ p_resid_trans, data = simulation_result)
p1 <- ggplot(simulation_result, aes(x = p_resid, y = p_lrt, col = label)) + 
  geom_point(size = 0.9) +
  theme_bw() +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("p (perm. test w/ residual statistic)") +
  ylab("p (GLM-based LRT)") +
  scale_color_manual(values = cols) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  geom_abline(intercept = fit_1$coefficients[1], slope = fit_1$coefficients[2],
              col = "darkorange", linetype = "dashed")
legend <- get_legend(p1)
p1 <- p1 + my_theme_no_legend

fit_2 <- lm(p_lrt_trans ~ p_score_trans, data = simulation_result)
p2 <- ggplot(simulation_result, aes(x = p_score, y = p_lrt, col = null_true)) + 
  geom_point(size = 0.9) +
  theme_bw() +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("p (perm. test w/ score statistic)") +
  ylab("p (GLM-based LRT)") +
  scale_color_manual(values = cols) + my_theme_no_legend +
  geom_abline(intercept = fit_2$coefficients[1], slope = fit_2$coefficients[2],
              col = "darkorange", linetype = "dashed")
p_final <- plot_grid(plot_grid(p1, p2, nrow = 1, labels = "auto"), legend, nrow = 2, rel_heights = c(0.9, 0.1))
to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_s8_s9/fig_s8.png")
ggsave(plot = p_final, filename = to_save_fp, device = "png",
       scale = 0.9, width = 6.5, height = 3.5, dpi = 330)

###########################
# PLOT 2: REAL DATA RESULTS
###########################
# function for analyzing discovery results
analyze_discovery_results <- function(dataset) {
  # set file paths
  full_base_dir <- paste0(sceptre2_dir, "results/discovery_analyses/", dataset, "_full_stat")
  resid_base_dir <- paste0(sceptre2_dir, "results/discovery_analyses/", dataset, "_resid_stat")
  
  # define several functions
  # define several functions
  join_result_dfs <- function(full_stat_df, resid_stat_df) {
    signif_present <- "significant" %in% colnames(full_stat_df)
    left_join(full_stat_df[,c("response_id", "grna_target", "p_value", if (signif_present) "significant" else NULL)],
              resid_stat_df[,c("response_id", "grna_target", "p_value", if (signif_present) "significant" else NULL)],
              by = c("response_id", "grna_target"),
              suffix = c("_full_stat", "_resid_stat"))
  }
  
  obtain_combined_df <- function(analysis_type, full_base_dir, resid_base_dir) {
    full_stat_df <- readRDS(paste0(full_base_dir, "/results_", analysis_type, ".rds")) |> as.data.frame()
    resid_stat_df <- readRDS(paste0(resid_base_dir, "/results_", analysis_type, ".rds")) |> as.data.frame()
    join_result_dfs(full_stat_df, resid_stat_df)
  }
  
  clip_result_df <- function(combined_discovery_analysis_df) {
    combined_discovery_analysis_df |>
      na.omit() |>
      mutate(p_value_full_stat = ifelse(p_value_full_stat < 1e-40, 1e-40, p_value_full_stat),
             p_value_resid_stat = ifelse(p_value_resid_stat < 1e-40, 1e-40, p_value_resid_stat))
  }
  
  # load combined data frames
  combined_calib_check_df <- obtain_combined_df("run_calibration_check", full_base_dir, resid_base_dir)
  combined_power_check_df <- obtain_combined_df("run_power_check", full_base_dir, resid_base_dir)
  combined_discovery_analysis_df <- obtain_combined_df("run_discovery_analysis", full_base_dir, resid_base_dir)

  # create plot
  p1 <- ggplot(clip_result_df(combined_calib_check_df),
               aes(y = p_value_full_stat, x = p_value_resid_stat)) +
    geom_point(size = 0.9, col = "indianred2", alpha = 0.7) +
    theme_bw() +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    geom_abline(slope = 1, intercept = 0) +
    ylab("p (score statistic)") +
    xlab("p (resid statistic)") +
    my_theme +
    ggtitle("Negative control pairs")
  
  p2 <- ggplot(combined_power_check_df |> na.omit(),
               aes(y = p_value_full_stat, x = p_value_resid_stat)) +
    geom_point(col = "mediumseagreen") +
    theme_bw() +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    geom_abline(slope = 1, intercept = 0) +
    ylab("p (score statistic)") +
    xlab("p (resid statistic)") +
    my_theme +
    ggtitle("Positive control pairs")
  
  p3 <- ggplot(data = combined_discovery_analysis_df |>
                 na.omit() |> filter(p_value_full_stat > 1e-249, p_value_resid_stat > 1e-249),
               aes(x = p_value_resid_stat, y = p_value_full_stat)) +
    geom_point(size = 0.9, col = "dodgerblue3", alpha = 0.7) +
    theme_bw() +
    scale_x_continuous(trans = revlog_trans(base = 10)) +
    scale_y_continuous(trans = revlog_trans(base = 10)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("p (resid statistic)") +
    ylab("p (score statistic)") +
    my_theme +
    ggtitle("Discovery pairs")
  
  p4 <- ggplot() + theme_minimal()
  p_final <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = "auto")
  
  # compute the number of discoveries (both discovery and calibration) as a function of FDR level 
  comute_n_rejections <- function(result_df) {
    fdr_levels <- c(0.1, 0.05, 0.01, 0.005)
    n_rejections_df <- result_df |>
      na.omit() |>
      select(p_value_full_stat, p_value_resid_stat) |>
      pivot_longer(cols = c("p_value_full_stat", "p_value_resid_stat"),
                   names_to = "method", values_to = "p_value") |>
      group_by(method) |>
      mutate(q_value = p.adjust(p_value, method = "BH")) |>
      reframe(n_reject = sapply(fdr_levels, function(fdr_level) sum(q_value < fdr_level)),
              fdr_level = fdr_levels)
    percent_change_df <- n_rejections_df |> 
      group_by(fdr_level) |>
      reframe(percent_change = 100 * (n_reject[method == "p_value_full_stat"]/n_reject[method == "p_value_resid_stat"] - 1))
    return(list(n_rejections_df = n_rejections_df, percent_change_df = percent_change_df))    
  }

  # obtain the running times
  resid_running_time <- readRDS(paste0(resid_base_dir, "/running_times.rds"))/(60^2)
  full_running_time <- readRDS(paste0(full_base_dir, "/running_times.rds"))/(60^2)
  running_time_matrix <- matrix(data = c(resid_running_time, full_running_time), nrow = 2)
  rownames(running_time_matrix) <- c("calibration_check", "discovery_analysis")
  colnames(running_time_matrix) <- c("resid", "full")
  
  # return all
  list(plot = p_final, running_time_matrix = running_time_matrix,
       n_rejections_calib = comute_n_rejections(combined_calib_check_df),
       n_rejections_discovery = comute_n_rejections(combined_discovery_analysis_df))
}

papalexi_res <- analyze_discovery_results("papalexi")
frangieh_res <- analyze_discovery_results("frangieh")

#to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
#                     "sceptre2-manuscript/R_scripts/figure_creation/fig_s8_s9/fig_s9.png")
#ggsave(filename = to_save_fp, plot = p_final, device = "png", width = 6, height = 5, dpi = 330)