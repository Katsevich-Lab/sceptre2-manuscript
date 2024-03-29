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
  "results/discovery_analyses/discovery_results_0423_processed.rds"
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
              mutate(Method = "Truth", odds_ratio = NA, pvalue = NA))

p_rejections <- results |>
  ggplot(aes(x = num_rejections, y = TF, fill = Method)) +
  geom_col(position = "dodge", color = "black") +
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 1000)) +
  scale_fill_manual(values = my_values) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  labs(x = "Number of rejections",
       y = "TF")

p_odds_ratios <- results |>
  ggplot(aes(x = odds_ratio, y = TF, fill = Method)) +
  geom_col(position = "dodge", color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4)) +
  scale_fill_manual(values = my_values) +
  theme(legend.position = "none") +
  labs(x = "Odds ratio",
       y = NULL)

p_pvalues <- results |>
  ggplot(aes(x = pvalue, y = TF, fill = Method)) +
  geom_col(position = "dodge", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(trans = katlabutils::revlog_trans(),
                     breaks = c(1e-20, 1e-40, 1e-60)) +
  scale_fill_manual(values = my_values) +
  theme(legend.position = "") +
  labs(x = "p-value",
       y = NULL)

if (FALSE) { # no need for fdp/power
  p_fdp <- results |>
    ggplot(aes(x = fdp, y = TF, fill = Method)) +
    geom_col(position = "dodge", color = "black") +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    scale_x_continuous(limits = c(0,1)) +
    scale_fill_manual(values = my_values) +
    theme(legend.position = "none") +
    labs(x = "FDP",
         y = NULL)
  
  p_power <- results |>
    ggplot(aes(x = power, y = TF, fill = Method)) +
    geom_col(position = "dodge", color = "black") +
    scale_fill_manual(values = my_values) +
    theme(legend.position = "") +
    labs(x = "Power",
         y = NULL)  
}

plot_grid(plot_grid(p_rejections + theme(legend.position = "none"),
                    p_fdp,
                    p_power,
                    p_odds_ratios,
                    p_pvalues,
                    nrow = 1),
          get_legend(p_rejections),
          ncol = 1,
          axis = "tlrb",
          rel_heights = c(1, 0.2))

dev.off()