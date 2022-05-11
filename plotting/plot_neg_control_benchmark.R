calibration_result <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                             "results/undercover_gRNA_check_results.rds") |> readRDS() |>
  dplyr::mutate(dataset = factor(x = dataset, levels = c("schraivogel_tap", "papalexi_gene", "schraivogel_perturb"), labels = c("Schraivogel TAP", "Papalexi Gene", "Schraivogel Perturb")),
                method = factor(x = method, levels = c("schraivogel_method", "seurat_de"),  labels = c("Schraivogel Method", "Seurat DE")))
code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/")

library(katlabutils)
library(ggplot2)
library(cowplot)

# check number of pairs analyzed by each method on each dataset
calibration_plot <- calibration_result |>
  ggplot(mapping = aes(y = p_value)) +
  facet_grid(method ~ dataset, scales = "free_x") +
  stat_qq_points(size = 0.7, ymin = 1e-10) +
  stat_qq_band() +
  geom_abline() +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  labs(x = "Expected quantile", y = "Observed quantile") +
  theme_bw()

ggsave(filename = paste0(code_dir, "figures/undercover_gRNA.png"), plot = calibration_plot, device = "png", scale = 1.2, width = 6, height = 4, dpi = 320)

# next, make a plot of RAM and CPU usage
to_plot <- calibration_result |>
  dplyr::select(undercover_gRNA, dataset, method, clock_time, max_ram) |>
  dplyr::distinct() |>
  dplyr::group_by(dataset, method) |>
  dplyr::summarize(m_clock_time = mean(clock_time)/60, m_max_ram = mean(max_ram)) |>
  dplyr::ungroup() |>
  tidyr::pivot_longer(cols = c("m_clock_time", "m_max_ram"),
                      names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = factor(x = metric, levels = c("m_clock_time", "m_max_ram"), labels = c("Time (m)", "RAM (GB)")),
                dataset = factor(x = dataset, levels = c("schraivogel_tap", "papalexi_gene", "schraivogel_perturb"),
                                 labels = c("Schraivogel TAP", "Papalexi gene", "Schraivogel Perturb"))) |>
  dplyr::mutate(value = ifelse(value < 1, 1.05, value))

out_p2 <- ggplot(to_plot, mapping = aes(x = dataset, y = value, fill = method)) + geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(.~metric, scales = "free_y") + ylab("") + xlab("Dataset") + scale_x_discrete(guide = guide_axis(n.dodge = 2)) + theme_bw() + scale_y_continuous(trans='log2')

ggsave(filename = paste0(code_dir, "figures/undercover_gRNA_compute.png"), plot = out_p2, device = "png", scale = 1.2, width = 6, height = 3, dpi = 320)
