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
# load simulation results, do basic preprocessing
simulation_result_fps <- paste0(sceptre2_dir, "results/extra_analyses/resid_vs_score_sim/") |> list.files(full.names = TRUE)
result_df <- lapply(X = seq_along(simulation_result_fps), FUN = function(i) {
  fp <- simulation_result_fps[i]
  readRDS(fp) |> mutate(proc_id = i)
}) |> data.table::rbindlist() |>
  mutate(label = ifelse(null_true, "Null true", "Alternative true"))

# apply BH, and obtain the mean false discovery proportion (fdr-hat) across datasets
fdr_level <- 0.1
cols <- c("deepskyblue4", "firebrick2")
n_nonnull <- result_df |>
  filter(proc_id == 1) |>
  summarize(sum(label == "Alternative true")) |>
  pull()

stat_summary_df <- result_df |> select(p_resid, p_score, p_lrt, null_true, proc_id) |> 
  pivot_longer(cols = c("p_resid", "p_score", "p_lrt"),
               names_to = "method", values_to = "p_val") |>
  group_by(method, proc_id) |>
  mutate(p_adj = p.adjust(p_val, method = "BH"), signif = p_adj < fdr_level) |>
  filter(signif) |>
  summarize(n_total_discoveries = dplyr::n(),
            n_false_discoveries = sum(null_true),
            fdp = n_false_discoveries/n_total_discoveries) |>
  summarize(mean_fdp = mean(fdp),
            mean_n_discoveries = mean(n_total_discoveries))
time_summary_df <- result_df |> select(resid_time, score_time, lrt_time) |>
  pivot_longer(cols = c("resid_time", "score_time", "lrt_time"),
               names_to = "method", values_to = "time") |>
  group_by(method) |>
  summarize(m_time = mean(time), 
            upper_ci = m_time + 1.96 * sd(time)/sqrt(dplyr::n()),
            lower_ci = m_time - 1.96 * sd(time)/sqrt(dplyr::n()))

# create plot from one run
curr_proc_id <- 7
simulation_result <- result_df |>
  filter(proc_id == curr_proc_id) |>
  mutate(p_resid_trans = -log(p_resid, 10),
         p_score_trans = -log(p_score, 10),
         p_lrt_trans = -log(p_lrt, 10))

my_breaks <- 10^(seq(from = 0, to = -6, by = -2))
fit_1 <- lm(p_lrt_trans ~ p_resid_trans, data = simulation_result)
p1 <- ggplot(simulation_result, aes(x = p_resid, y = p_lrt, col = label)) + 
  geom_point(size = 0.9) +
  theme_bw() +
  scale_x_continuous(trans = revlog_trans(base = 10), breaks = my_breaks) +
  scale_y_continuous(trans = revlog_trans(base = 10), breaks = my_breaks) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("p (perm. residual statistic)") +
  ylab("p (GLM-based LRT)") +
  scale_color_manual(values = cols) + my_theme_no_legend +
  geom_abline(intercept = fit_1$coefficients[1], slope = fit_1$coefficients[2],
              col = "darkorange", linetype = "dashed")

fit_2 <- lm(p_lrt_trans ~ p_score_trans, data = simulation_result)
p2 <- ggplot(simulation_result, aes(x = p_score, y = p_lrt, col = label)) + 
  geom_point(size = 0.9) +
  theme_bw() +
  scale_x_continuous(trans = revlog_trans(base = 10), breaks = my_breaks) +
  scale_y_continuous(trans = revlog_trans(base = 10), breaks = my_breaks) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("p (perm. score statistic)") +
  ylab("p (GLM-based LRT)") +
  scale_color_manual(values = cols) + my_theme +
  geom_abline(intercept = fit_2$coefficients[1], slope = fit_2$coefficients[2],
              col = "darkorange", linetype = "dashed") +
  theme(strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.1),
        legend.key.size = unit(0.35, 'cm'),
        legend.margin = margin(t = -0.5, unit = 'cm'),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))

fit_3 <- lm(p_score_trans ~ p_resid_trans, data = simulation_result)
p3 <- ggplot(simulation_result, aes(x = p_resid, y = p_score, col = label)) + 
  geom_point(size = 0.9) +
  theme_bw() +
  scale_x_continuous(trans = revlog_trans(base = 10), breaks = my_breaks) +
  scale_y_continuous(trans = revlog_trans(base = 10), breaks = my_breaks) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("p (perm. residual statistic)") +
  ylab("p (perm. score statistic)") +
  scale_color_manual(values = cols) + my_theme_no_legend +
  geom_abline(intercept = fit_3$coefficients[1], slope = fit_3$coefficients[2],
              col = "darkorange", linetype = "dashed")
p4 <- ggplot() + theme_minimal()
p_final <- plot_grid(plot_grid(p1, p2, p3, p4, nrow = 2, labels = "auto"))
to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_s8_s9/fig_s8.png")

ggsave(plot = p_final, filename = to_save_fp, device = "png",
       scale = 1.0, width = 6.5, height = 5.25, dpi = 330)