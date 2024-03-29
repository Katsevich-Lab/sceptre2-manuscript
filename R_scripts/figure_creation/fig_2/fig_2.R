# a) Four seurat DE resampling distributions on IFN gamma data, with N(0,1) density superimposed
# b) y-axis: p-value ratio, x-axis, KS fit quality, point color: sample size, point: Frangieh IFN-gamma pair
# c) Frangieh IFN gamma: Seurat w/ QC, Seurat w/o QC, NB reg w/ QC, NB reg w/o QC (legend in Figure)
# d) confounding on Papalexi data
# e) Papalexi: NB w/ cov, NB w/o cov (legend in plot)
# f) model misspecification on Frangieh and Papalexi (e.g., overlayed histograms of goodness of fit p-values)

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
                                 "undercover_grna_analysis/undercover_result_grp_1_0523_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF)
resampling_res <- readRDS(paste0(result_dir, "resampling_distributions/seurat_resampling_at_scale_processed.rds")) |>
  mutate(p_rat = p_emp/p_value)
fisher_exact_p <- readRDS(paste0(result_dir, "extra_analyses/papalexi_grna_confounding_tests.rds"))

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
  lab <- paste0("Pair ", i, " (KS stat = ", round(ks_stat, 3), ")")
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
        legend.position = c(0.29, 0.78),
        legend.key.size = unit(0.35, 'cm'),
        legend.margin = margin(t = -0.5, unit = 'cm'),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))
        ) +
  ylab("Density") +
  scale_y_continuous(expand = expansion(mult = c(0.0, .01))) +
  geom_line(aes(x = z_grid, y = density, col = "N(0,1) density"),
            data = density_df, linewidth = 0.7) +
  ggtitle("Permutation distribution of Wilcox. statistic") +
  xlab("Permuted Wilcox. statistic") +
  scale_color_manual(values = c("N(0,1) density" = "purple3")) +
  scale_fill_manual(values = c("Permutation distribution" = "lightgrey"))

##########
# PANNEL b
##########
pairs_to_annotate <- right_join(x = resampling_res,
                                y = pairs,
                                by = c("undercover_grna", "response_id")) |> arrange(pair_id)

p_b <- ggplot(data = resampling_res |> filter(p_rat < 10, p_rat > 1e-3, n_nonzero_treatment >= 1),
       mapping = aes(x = ks_stat, y = p_rat, col = log(n_nonzero_treatment + 1))) + 
  geom_point(alpha = 0.7, size = 0.8) +
  scale_y_log10() +
  scale_x_log10() +
  ylab(expression(italic(p)[ratio]*" = "*italic(p)[exact]/italic(p)[asymptotic])) +
  xlab("KS statistic") +
  geom_hline(yintercept = 1) +
  my_theme +
  theme(legend.position = c(0.1, 0.67),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank(),
        legend.margin=margin(t = -0.3, unit='cm')) +
  ggtitle("Inflation of Wilcox. p-values") +
  annotate("text", x = 0.007, y = 9.5, label = "Log(Effective sample size + 1)", size = 3) +
  # annotate pair 1
  geom_segment(aes(x = pairs_to_annotate[1,"ks_stat"],
                   xend = pairs_to_annotate[1,"ks_stat"],
                   yend = pairs_to_annotate[1,"p_rat"] + 0.13,
                   y = pairs_to_annotate[1, "p_rat"] + 2),
               col = "black",
               arrow = arrow(length = unit(0.03, "npc")),
               linewidth = 0.45) +
  annotate(geom = "label",
           x = pairs_to_annotate[1,"ks_stat"],
           y = pairs_to_annotate[1, "p_rat"] + 2,
           label = "Pair 1",
           size = 3) +
  # annotate pair 2
  geom_segment(aes(x = pairs_to_annotate[2,"ks_stat"] - 0.12,
                   xend = pairs_to_annotate[2,"ks_stat"] - 0.015,
                   yend = pairs_to_annotate[2,"p_rat"],
                   y = pairs_to_annotate[2, "p_rat"]),
               col = "black",
               arrow = arrow(length = unit(0.03, "npc")),
               linewidth = 0.45) +
  annotate(geom = "label",
           x = pairs_to_annotate[2,"ks_stat"] - 0.14,
           y = pairs_to_annotate[2, "p_rat"],
           label = "Pair 2",
           size = 3) +
  scale_color_gradient(low = "mediumpurple3",
                       high = "navy")
##############
# New pannel c
##############
pannel_c_cols <- c("darkorchid4", "dodgerblue4", "dodgerblue3", "dodgerblue1", "deepskyblue") 
frangieh_binned_undercover <- undercover_res |>
  filter(method == "seurat_de",
         dataset == "frangieh_ifn_gamma_gene") |>
  dplyr::mutate(n_nonzero_trt_bin = bin(n_nonzero_treatment, 5))
# compute the minimum number of pairs per bin
n_to_sample <- frangieh_binned_undercover |>
  group_by(n_nonzero_trt_bin) |>
  summarize(count = n()) |>
  pull(count) |> min()
# downsample
set.seed(1)
to_plot_c <- frangieh_binned_undercover |>
  group_by(n_nonzero_trt_bin) |>
  sample_n(n_to_sample)
p_c <- ggplot(data = to_plot_c |> dplyr::arrange(n_nonzero_trt_bin),
              mapping = aes(y = p_value, col = n_nonzero_trt_bin)) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10)) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  my_theme +
  annotate("text", x = 9e-2, y = 1e-7, label = "Effective sample size", size = 3) +
  theme(legend.position=c(0.12, 0.55),
        legend.margin=margin(t = -0.5, r = -0.5, unit = 'cm'),
        legend.title = element_blank()) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_c_cols, name = "N treatment cells with nonzero expression") +
  ggtitle("Wilcox. calibration on Frangieh IFN-\u03B3")

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
                               top = ggpubr::text_grob("Confounding (Papalexi gene modality)",
                                                       size = 11, hjust = 0.395))
###########
# PANNEL E
###########
my_labels <- c("NB regression (w/ covariates)", "NB regression (no covariates)")
undercover_res_sub <- undercover_res |>
  filter(method %in% c("nb_regression_w_covariates", "nb_regression_no_covariates"),
         dataset == "papalexi_eccite_screen_gene",
         n_nonzero_treatment >= 10,
         n_nonzero_control >= 10)
p_e <- undercover_res_sub |>
  ggplot(aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10)) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  my_theme +
  theme(legend.title= element_blank(),
        legend.position = c(0.38, 0.75),
        legend.margin=margin(t = -0.5, unit = 'cm')) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = my_cols[names(my_cols) %in% my_labels]) +
  ggtitle("Papalexi (gene) negative control pairs")

##########
# PANNEL F
##########
pannel_f_cols <- RColorBrewer::brewer.pal(6, "YlOrRd")[seq(6,2)]
frangieh_binned_undercover <- undercover_res |> 
  filter(method == "nb_regression_w_covariates",
         dataset == "frangieh_ifn_gamma_gene") |>
  dplyr::mutate(n_nonzero_trt_bin = bin(n_nonzero_treatment, 5))
# compute the minimum number of pairs per bin
n_to_sample <- frangieh_binned_undercover |>
  group_by(n_nonzero_trt_bin) |>
  summarize(count = n()) |>
  pull(count) |> min()
# downsample
set.seed(1)
to_plot_f <- frangieh_binned_undercover |>
  group_by(n_nonzero_trt_bin) |>
  sample_n(n_to_sample)
p_f <- ggplot(data = to_plot_f |> dplyr::arrange(n_nonzero_trt_bin),
              mapping = aes(y = p_value, col = n_nonzero_trt_bin)) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10)) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  my_theme +
  annotate("text", x = 9e-2, y = 1e-7, label = "Effective sample size", size = 3) +
  theme(legend.position=c(0.12, 0.55),
        legend.margin=margin(t = -0.5, r = -0.5, unit = 'cm'),
        legend.title = element_blank()) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_f_cols, name = "N treatment cells with nonzero expression") +
  ggtitle("NB reg. calibration on Frangieh IFN-\u03B3")

############
# CREATE FIG
############
fig <- cowplot::plot_grid(p_a, p_b,
                          p_c, p_d,
                          p_e, p_f, ncol = 2, labels = "auto", align = "vh", axis = "l")
to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_2/fig_2.png")
ggsave(filename = to_save_fp, plot = fig, device = "png",
       scale = 1.1, width = 6.5, height = 7, dpi = 330)
