calibration_result <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                             "undercover_gRNA_check_results.rds") |> readRDS()
code_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/")

library(utilities)
library(ggplot2)
library(cowplot)

# first, seurat method
seurat_p <- calibration_result |> dplyr::filter(method == "seurat_de") |>
  ggplot(mapping = aes(y = p_value)) +
  stat_qq_points(size = 0.7) +
  stat_qq_band() +
  geom_abline() +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  facet_wrap(~dataset, scales = "free") +
  labs(x = "Expected quantile", y = "Observed quantile") +
  theme_bw() +
  ggtitle("Seurat DE method")

# next, schraivogel method
schraivogel_p <- calibration_result |> dplyr::filter(method == "schraivogel_method") |>
  ggplot(mapping = aes(y = p_value)) +
  stat_qq_points(size = 0.7) +
  stat_qq_band() +
  geom_abline() +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  facet_wrap(~dataset, scales = "free") +
  labs(x = "Expected quantile", y = "Observed quantile") +
  theme_bw() +
  ggtitle("Schraivogel MAST method")

out_p <- cowplot::plot_grid(seurat_p, schraivogel_p, labels = c("a", "b"), nrow = 2)
ggsave(filename = paste0(code_dir, "figures/undercover_gRNA.png"), plot = out_p, device = "png", scale = 1.2, width = 6, height = 4, dpi = 320)
