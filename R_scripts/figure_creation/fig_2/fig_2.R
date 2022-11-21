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

shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
undercover_res <- readRDS(paste0(result_dir, "undercover_grna_analysis/undercover_result_grp_1_processed.rds"))

##########
# PANNEL a
##########
pairs <- data.frame(undercover_grna = c("NO-SITE-706", "NO-SITE-706", "NO-SITE-706", "NO-SITE-706"),
                    response_id = c("A1BG", "A1BG-AS1", "A4GALT", "AAAS"))

B <- 1000
pannel_a_list <- lapply(X = seq(1, nrow(pairs)), FUN = function(i) {
  args_to_pass <- lowmoi::get_sceptre_function_args_for_pair(response_id = pairs$response_id[i],
                                                             undercover_grna = pairs$undercover_grna[i],
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
  interval <- c(-3.75, 3.75)
  z_grid <- seq(interval[1], interval[2], length.out = 1000)
  density_df <- data.frame(z_grid = z_grid,
                           density = stats::dnorm(z_grid)) |>
    mutate(pair = paste0("pair_", i),
           ks_fit = paste0("Goodness of fit = ", round(res$ks_stat, 3)))
  histogram_df <- data.frame(z_null = z_null) |>
    mutate(pair = paste0("pair_", i),
           ks_fit = paste0("Goodness of fit = ", round(res$ks_stat, 3)))
  
  return(list(density_df = density_df, histogram_df = histogram_df))
})

density_df <- lapply(pannel_a_list, function(l) l$density_df) |>
  bind_rows() |>
  mutate(ks_fit = factor(ks_fit))
histogram_df <- lapply(pannel_a_list, function(l) l$histogram_df) |>
  bind_rows() |>
  mutate(ks_fit = factor(ks_fit))

p_a <- ggplot() + geom_histogram(aes(x = z_null, y = after_stat(density)),
                                 data = histogram_df,
                                 boundary = 0, color = "black",
                                 fill = "lightgray", bins = 15) +
  facet_wrap(ks_fit ~ ., nrow = 2, scales = "free_y") +
  my_theme +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank()) +
  ylab("") +
  scale_y_continuous(expand = expansion(mult = c(0.0, .01))) +
  geom_line(aes(x = z_grid, y = density),
            data = density_df, col = "darkred", lwd = 0.7) +
  ggtitle("Empirical null distribution of MW statistic")

##########
# PANNEL b
##########

##########
# PANNEL c
##########
# filter for seurat and NB reg (with cov) on Frangieh IFN gamma; add column for pass stringent QC
my_labels <- c("Seurat De (w/ strict QC)", "Seurat De (w/o strict QC)", "NB Reg (w/ strict QC)", "NB Reg (w/o strict QC)")
undercover_res_sub <- undercover_res |>
  filter(method %in% c("seurat_de", "nb_regression_w_covariates"),
         dataset == "frangieh_ifn_gamma_gene") |>
  mutate(pass_stringent_qc = (n_nonzero_treatment >= 30 & n_nonzero_control >= 30),
         Method = fct_recode(Method, "NB Reg"= "Nb Regression W Covariates")) |>
  mutate(Method = paste0(Method, ifelse(pass_stringent_qc, " (w/ strict QC)", " (w/o strict QC)"))) |>
  mutate(Method = fct_relevel(Method,
                              my_labels))
# compute the minimum number across method-QC status pairs
n_to_sample <- undercover_res_sub |>
  group_by(Method) |>
  summarize(count = n()) |>
  pull(count) |> min()
# downsample
to_plot_c <- undercover_res_sub |>
   group_by(Method) |>
   sample_n(n_to_sample)
p_c <- ggplot(data = to_plot_c, mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.55) +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10)) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  my_theme +
  theme(legend.title= element_blank(),
        legend.position = c(0.7, 0.12),
        legend.margin=margin(t = -0.5, unit='cm')) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.1,
    default.unit = "inch")) +
  scale_color_manual(values = my_cols[names(my_cols) %in% my_labels]) +
  ggtitle("Frangieh IFN-\u03B3 negative control pairs")

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
fisher_exact_p <- sapply(unique_nt_grnas, function(unique_nt_grna) {
  grna_binary_vect <- as.integer(unique_nt_grna == nt_grna_assignments)
  cont_table <- as.matrix(table(biorep_vect_nt, grna_binary_vect))
  fisher.test(x = cont_table)$p.value
})
my_nt_grna <- names(which.min(fisher_exact_p))
grna_binary_vect <- as.integer(my_nt_grna == nt_grna_assignments)
cont_table <- as.matrix(table(biorep_vect_nt, grna_binary_vect))
prop_table <- cont_table |>
  prop.table(margin = 1) |>
  as.data.frame() |>
  filter(grna_binary_vect == 1) |>
  select(-grna_binary_vect, freq = Freq, bio_rep = biorep_vect_nt) |>
  mutate(bio_rep = fct_recode(bio_rep, "Rep 1" = "rep_1", "Rep 2" = "rep_2", "Rep 3" = "rep_3"))

# carry out a similar analysis for relative gene expression
gene_exp_mat <- as.matrix(response_odm[[,nt_cells]])
rownames(gene_exp_mat) <- get_feature_ids(response_odm)
cell_cov <- (response_odm |> get_cell_covariates())[nt_cells,]
gene_ids <- response_odm |>
  get_feature_covariates() |>
  arrange(desc(mean_expression)) |>
  slice(1:50) |>
  row.names()
full_formula <- formula(expressions ~ bio_rep + offset(log(n_umis)))
reduced_formula <-  formula(expressions ~ offset(log(n_umis)))
lrt_p <- sapply(X = gene_ids, function(gene_id) {
  print(paste0("Fitting model for ", gene_id))
  expressions <- gene_exp_mat[gene_id,]
  curr_cell_cov <- mutate(cell_cov, expressions = expressions)
  # estimate the size
  est_size <- lowmoi:::estimate_size(df = curr_cell_cov,
                                     formula = full_formula)
  # fit the NB regression models
  full_nb_reg <- glm(formula = full_formula,
                     family = MASS::negative.binomial(est_size),
                     data = curr_cell_cov)
  reduced_nb_reg <- glm(formula = reduced_formula,
                        family = MASS::negative.binomial(est_size),
                        data = curr_cell_cov)
  fit <- anova(reduced_nb_reg, full_nb_reg, test = "LRT")
  fit$`Pr(>Chi)`[2]
}) # All LRTs are highly significant, indicating a strong association between (relative) gene expression and biological replicate. I will (somewhat arbitrarily) choose the second most highly expressed gene to plot.
gene_to_plot <- gene_ids[2] 
rel_expression_df <- data.frame(rel_expression = 1000 * log(gene_exp_mat[gene_to_plot,]/cell_cov[,"n_umis"] + 1),
                                bio_rep = cell_cov[,"bio_rep"]) |>
  mutate(bio_rep = fct_recode(bio_rep, "Rep 1" = "rep_1", "Rep 2" = "rep_2", "Rep 3" = "rep_3"))

p_d1 <- ggplot(data = prop_table,
               aes(x = bio_rep, y = freq)) +
  geom_bar(stat = "identity", fill = "lightgray", col = "black") +
  ylab("Frac. cells with perturbation") +
  xlab("Biological replicate") +
  my_theme +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0.01, 0.01))) +
  theme(plot.title = element_text(hjust=1))

p_d2 <- ggplot(data = rel_expression_df,
       aes(x = bio_rep, y = rel_expression)) +
  geom_violin(fill = "lightgrey", col = "black") +
  geom_boxplot(outlier.shape = NA, coef = 0, fill = "lightgrey") +
  scale_y_continuous(limits = c(0, 50),
                     expand = expansion(mult = c(0.01, 0))) +
  ylab("Relative gene expression") +
  xlab("Biological replicate") +
  my_theme +
  theme(plot.title = element_text(hjust=1))

p_d <- gridExtra::grid.arrange(p_d1, p_d2, nrow = 1,
                               top = ggpubr::text_grob("Confounding (Papalexi gene modality)", size = 11))



############
# CREATE FIG
############
fig <- cowplot::plot_grid(p_a, NULL,
                          p_c, p_d,
                          NULL, NULL, ncol = 2, labels = "auto", align = "h", axis = "l")
to_save_fp <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "sceptre2-manuscript/R_scripts/figure_creation/fig_2/r_output.png")
ggsave(filename = to_save_fp, plot = fig, device = "png",
       scale = 1.1, width = 6.5, height = 7, dpi = 330)
