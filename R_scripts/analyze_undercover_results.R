# get fps/data
undercover_res_fp <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                            "results/undercover_grna_analysis/undercover_gRNA_check_results.rds")
undercover_res <- readRDS(undercover_res_fp)
fig_dir <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/figures/")

# load packages
library(ggplot2)
library(katlabutils)

# first, check the number of p-values for each dataset-method pair; confirm number of p-values coincides across methods for a given dataset
undercover_res |>
  dplyr::group_by(dataset, method) |>
  dplyr::summarize(count = dplyr::n()) |>
  dplyr::summarize(n_pvals_coincide = all(diff(count) == 0))

undercover_res |>
  dplyr::group_by(dataset, method) |>
  dplyr::summarize(count = sum(is.na(p_value)))

# for each dataset, plot the p-values of each method
res_proc <- undercover_res |>
  dplyr::group_by(dataset, method) |>
  dplyr::mutate(n_pairs = dplyr::n(),
                dataset_rename = stringr::str_to_title(gsub(pattern = "_",replacement = " ", x = dataset)),
                dataset_rename_w_pairs = paste0(dataset_rename, " (", n_pairs[1], " pairs)"),
                Method = stringr::str_to_title(gsub(pattern = "_",replacement = " ", x = method)))

# first, the p-value plot
p_undercover_stat <- ggplot(data = res_proc, mapping = aes(y = p_value, color = Method)) +
  geom_vline(xintercept = 0.01) +
  stat_qq_points(ymin = 1e-10) +
  facet_wrap(~dataset_rename_w_pairs, scales = "free", labeller = label_wrap_gen(35)) +
  geom_abline() +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  theme_bw() +
  labs(x = "Expected quantile", y = "Observed quantile") +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

comp_df <- res_proc |>
  dplyr::select(undercover_gRNA, Dataset = dataset_rename, Method, clock_time, max_ram) |>
  dplyr::distinct() |>
  dplyr::group_by(Dataset, Method) |>
  dplyr::summarize(m_clock_time = mean(clock_time)/60, m_max_ram = mean(max_ram)) |>
  dplyr::ungroup() |>
  tidyr::pivot_longer(cols = c("m_clock_time", "m_max_ram"),
                      names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = factor(x = metric, levels = c("m_clock_time", "m_max_ram"),
                                labels = c("Time (m)", "RAM (GB)"))) |>
  dplyr::mutate(value = ifelse(value < 1, 1.05, value))
my_fct_order <- comp_df |>
  dplyr::filter(Method == "Schraivogel Method", metric == "RAM (GB)") |>
  dplyr::arrange(value) |>
  dplyr::pull(Dataset)
comp_df <- comp_df |>
  dplyr::mutate(Dataset = factor(x = Dataset, levels = my_fct_order, labels = my_fct_order))

p_undercover_comp <- ggplot(data = comp_df, mapping = aes(x = Method, y = value, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", col = "black") +
  facet_wrap(.~metric, scales = "free_y") +
  ylab("") + xlab("Method") +
  theme_bw() +
  scale_y_continuous(trans='log2')

ggsave(filename = paste0(fig_dir, "undercov_grna_stat.png"),
       plot = p_undercover_stat, device = "png", scale = 1, width = 11, height = 6, dpi = 330)
ggsave(filename = paste0(fig_dir, "undercov_grna_comp.pdf"),
       plot = p_undercover_comp, device = "pdf", scale = 0.8, width = 11, height = 4, dpi = 330)
