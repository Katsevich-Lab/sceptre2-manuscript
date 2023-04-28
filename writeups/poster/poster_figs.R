setwd("/Users/timbarry/research_code/sceptre2-manuscript/writeups/poster")

# Load packages
library(tidyverse)
library(katlabutils)
library(cowplot)
library(ondisc)
conflict_prefer(name = "filter", winner = "dplyr")

# set colors (not loaded by)
bio_rep_cols <- c("R1" = "firebrick3", "R2" = "navy", "R3" = "purple3")
bio_rep_fills <-  c("R1" = "firebrick3", "R2" = "navy", "R3" = "purple3")
dataset_cols <- c("frangieh_co_culture_gene" = "dodgerblue3",
                  "frangieh_control_gene" = "purple3",
                  "frangieh_ifn_gamma_gene" = "navy",
                  "papalexi_eccite_screen_gene" = "firebrick3")

# load functions and data
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                            "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
undercover_res <- readRDS(paste0(result_dir,
                                 "undercover_grna_analysis/undercover_result_grp_1_0423_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF)
resampling_res <- readRDS(paste0(result_dir, "resampling_distributions/seurat_resampling_at_scale_processed.rds")) |>
  mutate(p_rat = p_emp/p_value)
fisher_exact_p <- readRDS(paste0(result_dir, "extra_analyses/papalexi_grna_confounding_tests.rds"))
nb_gof_tests <- readRDS(paste0(result_dir, "extra_analyses/goodness_of_fit_tests.rds")) |>
  mutate(theta = NULL)




##########
# PANNEL a
##########
pairs <- data.frame(undercover_grna = c("NO-SITE-836", "NO-SITE-599"),
                    response_id = c("DNAL1", "ZNF615"),
                    pair_id = paste0("Pair ", 1:2))
B <- 2500
pannel_a_list <- lapply(X = seq(1, nrow(pairs)), FUN = function(i) {
  response_id <- pairs$response_id[i]
  undercover_grna <- pairs$undercover_grna[i]
  args_to_pass <- lowmoi::get_sceptre_function_args_for_pair(response_id = response_id,
                                                             undercover_grna = undercover_grna,
                                                             dataset_name = "frangieh/co_culture/gene",
                                                             output_amount = 2, B = B)
  response_odm <- args_to_pass$mm_odm |> ondisc::get_modality("response")
  grna_odm <- args_to_pass$mm_odm |> ondisc::get_modality("grna")
  res <- lowmoi::mann_whitney_perm(response_odm = response_odm,
                                   grna_odm = grna_odm,
                                   response_grna_group_pairs = args_to_pass$response_grna_group_pairs,
                                   B = B,
                                   progress = TRUE,
                                   full_output = TRUE)
  z_null <- res |> select(matches("z_[0-9]+")) |> as.numeric()
  ks_stat <- resampling_res |>
    filter(response_id == !!response_id,
           undercover_grna == !!undercover_grna) |> pull(ks_stat)
  z_star <- resampling_res |>
    filter(response_id == !!response_id,
           undercover_grna == !!undercover_grna) |> pull(z_star)
  
  interval <- c(-3.75, 3.75)
  z_grid <- seq(interval[1], interval[2], length.out = 1000)
  lab <- paste0("Pair ", i)
  density_df <- data.frame(z_grid = z_grid,
                           density = stats::dnorm(z_grid)) |>
    mutate(ks_fit = lab, z_star = z_star)
  histogram_df <- data.frame(z_null = z_null) |>
    mutate(ks_fit = lab, z_star)
  
  return(list(density_df = density_df, histogram_df = histogram_df))
})

density_df <- lapply(pannel_a_list, function(l) l$density_df) |>
  bind_rows() |>
  mutate(ks_fit = factor(ks_fit))
histogram_df <- lapply(pannel_a_list, function(l) l$histogram_df) |>
  bind_rows() |>
  mutate(ks_fit = factor(ks_fit))

p_a <- ggplot() + geom_histogram(aes(x = z_null, y = after_stat(density),
                                     fill = "Permutation distribution"),
                                 data = histogram_df,
                                 boundary = 0, color = "black", bins = 20) +
  facet_wrap(ks_fit ~ ., nrow = 1) +
  my_theme +
  theme(strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.35, 0.78),
        legend.key.size = unit(0.35, 'cm'),
        legend.margin = margin(t = -0.5, unit = 'cm'),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))
  ) +
  ylab("Density") +
  scale_y_continuous(expand = expansion(mult = c(0.0, .01))) +
  geom_line(aes(x = z_grid, y = density, col = "N(0,1) density"),
            data = density_df, linewidth = 0.7) +
  ggtitle("Permutation dist. of Wilcox. statistic") +
  xlab("Permuted Wilcox. statistic") +
  scale_color_manual(values = c("N(0,1) density" = "purple3")) +
  scale_fill_manual(values = c("Permutation distribution" = "lightgrey"))

ggsave(filename = "wilcox_dist.png", plot =  p_a, device = "png", scale = 0.6, width = 5, height = 4, dpi = 330)

##########
# PANNEL d
##########
# load papalexi data
response_odm <- lowmoi::load_dataset_modality("papalexi/eccite_screen/gene")
grna_odm <- lowmoi::load_dataset_modality("papalexi/eccite_screen/grna_assignment")

# find a gRNA group that is strongly associated with bio_rep
grna_assignments <- lowmoi:::get_grna_assignments_via_max_op(grna_odm)
nt_cells <- grepl(pattern = "^NTg*", x = grna_assignments)
nt_grna_assignments <- grna_assignments[nt_cells]
unique_nt_grnas <- unique(nt_grna_assignments)
biorep_vect <- response_odm |>
  get_cell_covariates() |>
  pull(bio_rep)
biorep_vect_nt <- biorep_vect[nt_cells]
my_nt_grna <- names(which.min(fisher_exact_p))
grna_binary_vect <- as.integer(my_nt_grna == nt_grna_assignments)
cont_table <- as.matrix(table(biorep_vect_nt, grna_binary_vect))
prop_table <- cont_table |>
  prop.table(margin = 1) |>
  as.data.frame() |>
  filter(grna_binary_vect == 1) |>
  select(-grna_binary_vect, freq = Freq, bio_rep = biorep_vect_nt) |>
  mutate(bio_rep = fct_recode(bio_rep, "R1" = "rep_1", "R2" = "rep_2", "R3" = "rep_3"))

# carry out a similar analysis for relative gene expression
gene_exp_mat <- as.matrix(response_odm[[,nt_cells]])
rownames(gene_exp_mat) <- get_feature_ids(response_odm)
cell_cov <- (response_odm |> get_cell_covariates())[nt_cells,]
gene_to_plot <- response_odm |>
  get_feature_covariates() |>
  arrange(desc(mean_expression)) |>
  slice(2) |>
  row.names()
full_formula <- formula(expressions ~ bio_rep + offset(log(n_umis)))
reduced_formula <-  formula(expressions ~ offset(log(n_umis)))
rel_expression_df <- data.frame(rel_expression = 1000 * log(gene_exp_mat[gene_to_plot,]/cell_cov[,"n_umis"] + 1),
                                bio_rep = cell_cov[,"bio_rep"]) |>
  mutate(bio_rep = fct_recode(bio_rep, "R1" = "rep_1", "R2" = "rep_2", "R3" = "rep_3"))

p_d1 <- ggplot(data = prop_table,
               aes(x = bio_rep, y = freq, col = bio_rep, fill = bio_rep)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  ylab("Frac. cells with perturbation") +
  xlab("Biological rep.") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0.01, 0.01))) +
  my_theme + theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 12, unit = "pt"),
                   axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0, unit = "pt")),
                   legend.position = "none") +
  scale_color_manual(values = bio_rep_cols) +
  scale_fill_manual(values = bio_rep_fills)

p_d2 <- ggplot(data = rel_expression_df,
               aes(x = bio_rep, y = rel_expression, col = bio_rep, fill = bio_rep)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, coef = 0, fill = NA) +
  scale_y_continuous(limits = c(0, 50),
                     expand = expansion(mult = c(0.01, 0))) +
  ylab("Relative gene expression") +
  xlab("Biological rep.") +
  my_theme +
  theme(plot.title = element_text(hjust=1),
        legend.position = "none") +
  scale_color_manual(values = bio_rep_cols) +
  scale_fill_manual(values = bio_rep_fills)

p_d <- gridExtra::grid.arrange(p_d1, p_d2, nrow = 1,
                               top = ggpubr::text_grob("Confounding due to biological rep.",
                                                       size = 11, hjust = 0.395))
ggsave(filename = "confounding.png", plot =  p_d, device = "png", scale = 0.6, width = 5, height = 4.2, dpi = 330)


#############
# left pannel
#############
undercover_res <- readRDS(paste0(result_dir,
                                 "undercover_grna_analysis/undercover_result_grp_1_0423_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF)
alpha <- 0.1

make_figure_row <- function(dataset, print_legend) {
  my_methods <- c("KS test", "MAST", "MIMOSCA", "t-test", "Seurat-Wilcox", "Seurat-NB")
  my_values <- my_cols[names(my_cols) %in% my_methods]
  
  df_sub <- undercover_res |>
    filter(dataset == !!dataset,
           Method %in% my_methods) |>
    dplyr::arrange(Method)
  
  bonf_thresh <- alpha/(nrow(df_sub)/length(my_methods))
  
  p2 <- ggplot(data = df_sub, mapping = aes(y = p_value, col = Method)) +
    stat_qq_points(ymin = 1e-8, size = 0.55) +
    stat_qq_band() +
    scale_x_continuous(trans = revlog_trans(10)) +
    scale_y_continuous(trans = revlog_trans(10)) +
    labs(x = "Expected null p-value", y = "Observed p-value") +
    geom_abline(col = "black") + 
    geom_hline(yintercept = bonf_thresh, linetype = "dashed") +
    scale_color_manual(values = my_values)
  
  if (print_legend) {
    p2 <- p2 +
      my_theme +
      theme(legend.title = element_blank(),
            legend.position = c(0.78, 0.22),
            legend.text=element_text(size = 8.0),
            legend.margin=margin(t = 0, unit='cm'),
            legend.background = element_blank()) +
      guides(color = guide_legend(
        keywidth = 0.0,
        keyheight = 0.1,
        default.unit="inch",
        override.aes = list(size = 2.5)))
  } else {
    p2 <- p2 + my_theme_no_legend
  }
  
  bonf_reject_df <- compute_n_bonf_rejections(df_sub, alpha = alpha)
  
  p3 <- bonf_reject_df |>
    ggplot2::ggplot(ggplot2::aes(x = Method, y = n_reject, fill = Method)) +
    ggplot2::geom_col(col = "black") +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.2),
                       expand = c(0,0),
                       breaks = c(0, 1, 10, 100, 1000, 8000)) +
    ylab("False positives") +
    xlab("Method") + my_theme_no_legend +
    theme(axis.text.x = element_blank()) +
    scale_fill_manual(values = my_cols)
  
  p_row <- plot_grid(p2, p3, nrow = 1, labels = NULL,
                     align = "h", rel_widths = c(2,1))
  p_row
}

r1 <- make_figure_row(dataset = "papalexi_eccite_screen_gene", print_legend = TRUE)
ggsave(filename = "init_benchmark.png", plot = r1, device = "png", scale = 0.9, width = 5, height = 2.5, dpi = 330)


#################################################################
# Load packages and resolve namespace conflicts
#################################################################

# Load packages
library(tidyverse)
library(katlabutils)
library(ggpubr)
library(grid)
library(gridExtra)
library(gtable)

# Resolve namespace conflicts
conflicts_prefer(dplyr::filter)

#################################################################
# Create QQ plots
#################################################################

# read colors from my_cols
my_values <- my_cols[names(my_cols) %in% c("Seurat-Wilcox", "Seurat-NB", "SCEPTRE")]

# Frangieh QQ plot
qq_frangieh <- undercover_res |>
  mutate(Method = fct_relevel(Method, "Seurat-Wilcox", "SCEPTRE", after = Inf)) |>
  filter(dataset == "frangieh_ifn_gamma_gene",
         method %in% c("sceptre", "seurat_de", "seurat_de_nb")) |> 
  ggplot(mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-10, size = 0.85) +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10)) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  # ggtitle("Frangieh (IFN-\u03B3) neg. controls") +
  ggtitle("Frangieh (IFN-\u03B3) neg. controls") +
  scale_color_manual(values = my_values) + 
  my_theme +
  theme(legend.title = element_blank(),
        legend.position = c(0.3, 0.86),
        legend.text = element_text(size = 11),
        legend.margin = margin(t = 0, unit = 'cm')) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.2,
    default.unit = "inch",
    override.aes = list(size = 2.5)))

# Papalexi QQ plot
qq_papalexi <- undercover_res |>
  mutate(Method = fct_relevel(Method, "Seurat-Wilcox", "SCEPTRE", after = Inf)) |>
  filter(dataset == "papalexi_eccite_screen_gene",
         method %in% c("sceptre", "seurat_de", "seurat_de_nb")) |> 
  ggplot(mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-9, size = 0.85) +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10)) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  ggtitle("Papalexi (gene) neg. controls") +
  scale_color_manual(values = my_values) + 
  my_theme_no_legend

#################################################################
# Put the pieces together and save
#################################################################

# put the pieces together
final_plot <- ggarrange(qq_frangieh, qq_papalexi, nrow = 1)

# save the figure
ggsave(filename = "right_pannel_fig.png", 
       plot = final_plot,
       device = "png", 
       width = 7,
       height = 3.5,
       bg = "white",
       dpi = 330,
       scale = 0.85)
