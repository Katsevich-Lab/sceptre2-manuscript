library(tidyverse)
source(paste0(.get_config_path("LOCAL_CODE_DIR"),
              "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R"))
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

# load the dataset sparsity results; determine mean sparsity level of gene modality, protein modality, and TAP-seq
frac_n_nonzero_df <- readRDS(paste0(sceptre2_dir, "/results/dataset_sample_sizes/frac_entries_zero.rds"))
gene_modality_sparsity <- 
  c(frac_n_nonzero_df |> dplyr::filter(paper == "frangieh", dataset == "co_culture", modality == "gene") |> pull(prop_zero),
    frac_n_nonzero_df |> dplyr::filter(paper == "frangieh", dataset == "control", modality == "gene") |> pull(prop_zero),
    frac_n_nonzero_df |> dplyr::filter(paper == "frangieh", dataset == "ifn_gamma", modality == "gene") |> pull(prop_zero),
    frac_n_nonzero_df |> dplyr::filter(paper == "papalexi", dataset == "eccite_screen", modality == "gene") |> pull(prop_zero)) |> mean()
tap_modality_sparsity <-
  c(frac_n_nonzero_df |> dplyr::filter(paper == "schraivogel", dataset == "enhancer_screen_chr11") |> pull(prop_zero),
    frac_n_nonzero_df |> dplyr::filter(paper == "schraivogel", dataset == "enhancer_screen_chr8") |> pull(prop_zero)) |> 
  mean()
protein_modality_sparsity <- frac_n_nonzero_df |> dplyr::filter(paper == "papalexi", dataset == "eccite_screen", modality == "protein") |> pull(prop_zero)

# there are three parameters:
# 1. the number of nonzero treatment cells threshold (7)
# 2. the sparsity
# 3. the probability that a pair passes pairwise QC
# we hold fixed the number of nonzero treatment cells threshold
# we vary the sparisty over the grid; the first entry represents a standard single-cell experiment, the second a TAP-seq experiment, and the final an experiment with protein readout
# similarly, we vary the fraction of pairs surviving pairwise QC over the grid (0.5, ..., 0.95)
# we compute the number of cells per perturbation required such that x% of pairs survive pairwise QC

sparsity_levels <- c(gene_modality_sparsity, tap_modality_sparsity, protein_modality_sparsity) |> round(2)
fractions_pass <- seq(from = 0.5, to = 0.99, by = 0.01)
n_nonzero_trt_threshes <- c(7, 14)
df <- sapply(sparsity_levels, function(sparsity_level) {
    sapply(fractions_pass, function(fraction_pass) {
    sapply(n_nonzero_trt_threshes, function(n_nonzero_trt_thresh) {
      curr_n_cells <- 1L
      repeat {
        pass_qc <- qbinom(p = fraction_pass, size = curr_n_cells,
                          prob = 1 - sparsity_level, lower.tail = FALSE) >= n_nonzero_trt_thresh
        if (pass_qc) {
          break
        } else {
          curr_n_cells <- curr_n_cells + 1L
        }
      }
      data.frame(sparsity_level = sparsity_level, fraction_pass = fraction_pass,
                 n_cells = curr_n_cells, n_nonzero_trt_thresh = n_nonzero_trt_thresh)
    }, simplify = FALSE)
  }, simplify = FALSE)
}, simplify = FALSE) |> unlist(recursive = FALSE) |> unlist(recursive = FALSE) |>
  data.table::rbindlist() |>
  mutate(sparsity_level = factor(sparsity_level),
         n_nonzero_trt_thresh = paste0("N nonzero treat thresh = ", n_nonzero_trt_thresh) |> forcats::fct_relevel("N nonzero treat thresh = 7"))

# make plot
p <- ggplot(data = df, mapping = aes(x = fraction_pass * 100, y = n_cells, col = sparsity_level)) +
  geom_line(linewidth = 0.9) + theme_bw() + facet_grid(. ~ n_nonzero_trt_thresh) +
  scale_y_continuous(breaks = c(0, 10, 20, 40, 60, 80, 100, 120)) +
  xlab("Percent pairs passing QC") + ylab("N cells per perturbation") +
  labs(color = "Sparsity") + 
  theme(legend.key.size = unit(0.35, 'cm'),
        legend.margin = margin(t = -0.5, unit = 'cm'),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  ggtitle("N cells per perturbation vs.percent pairs passing QC") +
  scale_color_manual(values = c("dodgerblue3", "darkorchid3", "firebrick3"))

# how many needed cells needed to for 90% of pairs to pass QC for a standard gene expression assay?
df |> dplyr::filter(sparsity_level == 0.78, fraction_pass == 0.9) # 46 cells for standard QC threshold, and roughly twice as many -- 84 -- for the more stringent pairwise QC threshold.

# save the plot
to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_s13/fig_s13.png")
ggsave(filename = to_save_fp, plot = p, device = "png", width = 6, height = 3, dpi = 330)
