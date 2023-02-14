# Load packages
library(tidyverse)
library(katlabutils)
library(cowplot)
library(janitor)
library(kable)
library(kableExtra)

# Load scripts and results
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
undercover_res <- readRDS(paste0(result_dir,
                                 "undercover_grna_analysis/undercover_result_grp_1_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF)

pc_res <- readRDS(paste0(result_dir, "positive_control_analysis/pc_results_processed.rds"))
reject_thresh <- 1e-5


alpha <- 0.1
n_false_rejections <- undercover_res |>
  filter(!(Method %in% c(c("Nb Regression No Covariates", 
                           "Nb Regression W Covariates", 
                           "Sceptre No Covariates")))) |>
  group_by(dataset_rename, Method) |>
  summarize(reject = (p_value < alpha/n()),
                   Method = Method[1]) |>
  summarize(n_false_reject = sum(reject),
                   Method = Method[1]) |>
  ungroup()

#################################################################
# Create tables
#################################################################

### Type-I error ###

n_false_rejections_tab <- n_false_rejections |>
  pivot_wider(names_from = Method, values_from = n_false_reject) |>
  relocate("Sceptre", .after = "dataset_rename") |>
  rename(Dataset = dataset_rename)

n_false_rejections_tab |>
  bind_rows(
    n_false_rejections_tab |>
      summarise(across(-Dataset, mean)) |>
      mutate(Dataset = "Average")
  )
  
### Power ###

n_pc_reject_df <- pc_res |>
  filter(n_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_control >= N_NONZERO_CONTROL_CUTOFF) |>
  group_by(dataset_rename, Method) |>
  summarize(n_pc_reject = sum(p_value < reject_thresh),
            `PC pairs` = n(),
            Method = Method[1]) |>
  group_by(dataset_rename) |>
  # mutate(Method = forcats::fct_recode(Method,
  #                                     "SCEPTRE" = "Sceptre",
  #                                     "Liscovitch Meth." = "Liscovitch Method",
  #                                     "Schraivogel Meth." = "Schraivogel Method",
  #                                     "Weissman Meth." = "Weissman Method"),
  #        Method = forcats::fct_relevel(Method,
  #                                      "SCEPTRE",
  #                                      "Seurat De",
  #                                      "Weissman Meth.",
  #                                      "Mimosca",
  #                                      "Liscovitch Meth.",
  #                                      "Schraivogel Meth.")) |>
  left_join(n_false_rejections,
            by = c("dataset_rename", "Method")) |>
  filter(n_false_reject <= 10) |> 
  select(dataset_rename, Method, n_pc_reject, `PC pairs`) |>
  ungroup()

type_I_err_table_filename <- paste0(.get_config_path("LOCAL_CODE_DIR"), 
                           "sceptre2-manuscript/R_scripts/figure_creation/genes_figs/type_I_err_table.png")

table <- n_pc_reject_df |>
  pivot_wider(names_from = Method, values_from = n_pc_reject) |>
  relocate("Sceptre", .after = "dataset_rename") |>
  relocate(`PC pairs`, .after = `Schraivogel Method`) |>
  rename(Dataset = dataset_rename) 

table |>
  kable(format = "latex", 
        row.names = NA,
        booktabs = TRUE,
        linesep = "",
        digits = 2) |>
  save_kable(type_I_err_table_filename)

#################################################################
# Proportion of rejections
#################################################################

#prop_reject_df <- 
  pc_res |>
  filter(n_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_control >= N_NONZERO_CONTROL_CUTOFF) |>
  group_by(dataset, method) |>
  summarize(prop_pc_reject = mean(p_value < reject_thresh),
            Method = Method[1]) |>
  group_by(dataset) |>
  mutate(
         Method = forcats::fct_recode(Method,
                                      "SCEPTRE" = "Sceptre",
                                      "Liscovitch Meth." = "Liscovitch Method",
                                      "Schraivogel Meth." = "Schraivogel Method",
                                      "Weissman Meth." = "Weissman Method"),
         Method = forcats::fct_relevel(Method,
                                       "SCEPTRE",
                                       "Seurat De",
                                       "Weissman Meth.",
                                       "Mimosca",
                                       "Liscovitch Meth.",
                                       "Schraivogel Meth.")) |>
  left_join(n_false_rejections,
            by = c("dataset", "method")) |>
  filter(n_ntc_reject < 50) |>
  select(dataset, method, prop_pc_reject) |>
  ungroup() |>
  pivot_wider(names_from = "method", values_from = "prop_pc_reject")

#################################################################
# Try ROC curves
#################################################################

# join positive and negative control results
joined_res <- bind_rows(
  undercover_res |>
    select(p_value, dataset_rename, Method) |>
    mutate(type = "NTC"),
  pc_res |>
    select(p_value, dataset_rename, Method) |>
    mutate(type = "PC")
)

# create ROC curves
p <- joined_res |> 
  filter(dataset_rename != "Simulated Experiment 1 Gene", 
         !(Method %in% c("Nb Regression No Covariates", 
                         "Nb Regression W Covariates", 
                         "Sceptre No Covariates"))) |> 
  group_by(dataset_rename, Method) |> 
  arrange(p_value) |> 
  mutate(TPR = cumsum(type == "PC")/sum(type == "PC")) |> 
  arrange(desc(p_value)) |> 
  mutate(FPR = 1-cumsum(type == "NTC")/sum(type == "NTC")) |> 
  ggplot(aes(x = FPR, y = TPR, color = Method)) + 
  geom_line() + 
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
  geom_abline(linetype = "dashed") + 
  facet_wrap(~dataset_rename) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# save plot
roc_fig_filename <- paste0(.get_config_path("LOCAL_CODE_DIR"), 
                           "sceptre2-manuscript/R_scripts/figure_creation/genes_figs/roc.png")
ggsave(filename = roc_fig_filename, 
       plot = p, 
       device = "png",
       width = 6.5, 
       height = 6, 
       dpi = 330)

#################################################################
# Look at Papalexi positive controls in the context of Mixscape
#################################################################

# list of perturbations surviving Mixscape
ptrb_surviving_mixscape <- c("SMAD4", "STAT2", "JAK2", 
                             "STAT1", "IFNGR2", "IFNGR1", 
                             "IRF1", "BRD4", "SPI1", 
                             "CUL3", "MYC")

# SCEPTRE results
pc_res |>
  filter(n_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_control >= N_NONZERO_CONTROL_CUTOFF, 
         dataset == "papalexi_eccite_screen_gene", 
         Method %in% c("Sceptre")) |> 
  arrange(p_value) |> 
  filter(p_value < reject_thresh) |> 
  select(response_id, p_value) |> 
  mutate(survived_mixscape = response_id %in% ptrb_surviving_mixscape) |> 
  rename(gene = response_id, `sceptre p-val` = p_value)

# Seurat DE results
pc_res |>
  filter(n_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_control >= N_NONZERO_CONTROL_CUTOFF, 
         dataset == "papalexi_eccite_screen_gene", 
         Method %in% c("Seurat De")) |> 
  arrange(p_value) |> 
  filter(p_value < reject_thresh) |> 
  select(response_id, p_value) |> 
  mutate(survived_mixscape = response_id %in% ptrb_surviving_mixscape) |> 
  rename(gene = response_id, `Seurat DE p-val` = p_value)