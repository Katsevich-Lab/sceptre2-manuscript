<<<<<<< HEAD
pc_res <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
                 "results/positive_control_analysis/pc_results_sceptre_unfiltered_0523_processed.rds") |>
  readRDS() |> na.omit()
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), 
                            "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
library(ggplot2)
library(dplyr)
library(katlabutils)
library(cowplot)
source(shared_fig_script)
reject_thresh <- 1e-5 
conflicts_prefer(dplyr::filter)
N_NONZERO_TREATMENT_CUTOFF <- 7

make_p_val_vs_sample_size_plot <- function(pc_res_sub, ylab = TRUE, xlab = TRUE) {
  tit <- as.character(pc_res_sub$dataset_rename[1])
  p <- pc_res_sub |>
    mutate(reject = ifelse(p_value < reject_thresh, "Signif.", "Not signif.")) |>
    mutate(p_value = ifelse(p_value < 1e-7, 1e-7, p_value)) |>
    ggplot(mapping = aes(x = n_treatment, y = p_value, col = reject)) +
    annotate("rect", xmin = -Inf, xmax = N_NONZERO_TREATMENT_CUTOFF, ymin = 0, ymax = Inf, fill = "slategray1") +
    annotate("rect", xmin = N_NONZERO_TREATMENT_CUTOFF, xmax = Inf, ymin = 0, ymax = Inf, fill = "lightpink") +
    geom_hline(yintercept = reject_thresh, linetype = "dashed", col = "darkred") +
    geom_point(alpha = 0.85, size = 0.95) +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 2),
                       breaks = c(0, 1, 10, 100), expand = c(0.02, 0)) +
    scale_y_continuous(trans = revlog_trans(10)) +
    xlab(if (xlab) "Effective sample size" else "") +
    ylab("SCEPTRE p-value") +
    ggtitle(tit) +
    scale_color_manual(values = c("Signif." = "black", "Not signif." = "slategrey"),
                       breaks = c("Signif.", "Not signif.")) +
    my_theme_no_legend +
    (if (!ylab) theme(axis.title.y = element_blank()))
  return(p)
}

p_a <- make_p_val_vs_sample_size_plot(pc_res |> filter(dataset == "frangieh_co_culture_gene"),
                                     ylab = TRUE, xlab = FALSE)
p_b <- make_p_val_vs_sample_size_plot(pc_res |> filter(dataset == "frangieh_control_gene"),
                                      ylab = FALSE, xlab = TRUE)
p_c <- make_p_val_vs_sample_size_plot(pc_res |> filter(dataset == "frangieh_ifn_gamma_gene"),
                                      ylab = FALSE, xlab = FALSE)

fig_bottom <- plot_grid(p_a, p_b, p_c, nrow = 1, labels = c("a", "b", "c"))

fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
             "sceptre2-manuscript/R_scripts/figure_creation/fig_s8/fig_s8.png")
ggsave(filename = fp, plot = fig_bottom, device = "png", scale = 1, width = 7, height = 3, dpi = 320)
=======
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

# shared script
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                            "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
my_methods <- c("KS test", "MAST", "MIMOSCA", "t-test", "Seurat-Wilcox", "Seurat-NB", "SCEPTRE")
my_values <- c(my_cols[names(my_cols) %in% my_methods])

# analysis parameteres
N_NONZERO_CNTRL_THRESH <- 7
N_NONZERO_TRT_THRESH <- 7
q <- 0.1
#add file path of figure here
figure_file_path = ""

# read the results
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
method_results_fp <- paste0(
  sceptre2_dir,
  "results/discovery_analyses/discovery_results_0124_processed.rds"
)
TF_targets_fp <- paste0(
  sceptre2_dir,
  "results/discovery_analyses/TF_targets_papalexi_chipseq.rds"
)
method_results <- readRDS(method_results_fp) |> as_tibble()
TF_targets <- readRDS(TF_targets_fp) |> unique()

# join method results with TF targets, apply BH at level q
method_results <- method_results |>
  filter(dataset_rename == "Papalexi (Gene)",
         n_control >= N_NONZERO_CNTRL_THRESH,
         n_treatment >= N_NONZERO_TRT_THRESH) |>
  rename(gene = response_id, TF = grna_group) |>
  inner_join(TF_targets, by = c("gene", "TF")) |>
  group_by(TF, Method) |>
  #  group_by(Method) |>
  mutate(rejection = p.adjust(p_value, method = "fdr") <= q) |>
  ungroup()

# version of Fisher test that returns NA (rather than throwing an error)
# in degenerate cases
my.fisher.test <- function(x, y) {
  if(length(unique(x)) == 1 | length(unique(y)) == 1){
    list(p.value = NA, estimate = NA)
  } else{
    fisher.test(x, y, alternative = "greater")
  }
}

results <- method_results |>
  group_by(TF, Method) |>
  summarise(num_rejections = sum(rejection),
            odds_ratio = my.fisher.test(target, rejection)$estimate,
            pvalue = my.fisher.test(target, rejection)$p.value,
            fdp = sum(rejection & !target)/max(1, sum(rejection)),
            power = sum(rejection & target)/sum(target),
            .groups = "drop") |>
  bind_rows(TF_targets |>
              group_by(TF) |>
              summarise(num_rejections = sum(target)) |>
              mutate(Method = "ChIP-seq", odds_ratio = NA, pvalue = NA)) |>
  mutate(Method = forcats::fct_relevel(Method, "KS test",  "MIMOSCA", "t-test", "MAST", "Seurat-Wilcox", "SCEPTRE", "Seurat-NB"))

p_rejections <- results |>
  ggplot(aes(x = num_rejections, y = TF, fill = Method)) +
  geom_col(position = "dodge", color = "black") +
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 1000), expand = c(0, NA)) +
  scale_y_discrete(expand = c(0, NA)) +
  scale_fill_manual(values = c(my_values, "ChIP-seq" = "grey50")) +
  my_theme + theme(legend.position = "bottom") +
  labs(x = "N downstream genes",
       y = "Transcription factor")
l <- get_legend(p_rejections)
p_rejections <- p_rejections + my_theme_no_legend

p_odds_ratios <- results |>
  ggplot(aes(x = odds_ratio, y = TF, fill = Method)) +
  geom_col(position = "dodge", color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4), expand = c(0, NA)) +
  scale_y_discrete(expand = c(0, NA)) +
  scale_fill_manual(values = my_values) +
  my_theme_no_legend +
  labs(x = "Odds ratio", y = NULL)

p_pvalues <- results |>
  mutate(pvalue = ifelse(pvalue > 0.2, 0.2, pvalue)) |>
  ggplot(aes(x = pvalue, y = TF, fill = Method)) +
  geom_col(position = "dodge", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(trans = katlabutils::revlog_trans(),
                     breaks = c(1e-20, 1e-40, 1e-60), expand = c(0, NA)) +
  scale_y_discrete(expand = c(0, NA)) +
  scale_fill_manual(values = my_values) +
  my_theme_no_legend +
  labs(x = "p-value", y = NULL)

p_top <- plot_grid(p_rejections, p_odds_ratios, p_pvalues, nrow = 1)
p <- plot_grid(p_top, l, nrow = 2, rel_heights = c(0.85, 0.15))
to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/fig_s8/fig_s8.png")
ggsave(filename = to_save_fp, plot = p, device = "png",
       scale = 1.2, width = 6.0, height = 2.75, dpi = 330)
>>>>>>> 357d1b5 (update chip seq fig)
