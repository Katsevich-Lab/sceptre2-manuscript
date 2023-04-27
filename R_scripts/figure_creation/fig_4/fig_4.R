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
# Set analysis parameters
#################################################################

reject_thresh <- 1e-5   # threshold for rejection of positive controls
alpha <- 0.1            # target FWER level for negative controls
max_false_reject <- 50  # maximum false rejections to display power

#################################################################
# Load results
#################################################################

# source shared figure script
shared_fig_script <- paste0(
  .get_config_path("LOCAL_CODE_DIR"), 
  "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R"
)
source(shared_fig_script)

# directory with results
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")

# results of undercover analysis
undercover_res <- readRDS(paste0(
  result_dir,
  "undercover_grna_analysis/undercover_result_grp_1_0423_processed.rds"
)) |>
  filter(
    n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
    n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF
  )
# results of positive control analysis
pc_res <- readRDS(paste0(
  result_dir,
  "positive_control_analysis/pc_results_processed.rds"
)) |> 
  dplyr::filter(n_treatment >= N_NONZERO_TREATMENT_CUTOFF,
                n_control >= N_NONZERO_CONTROL_CUTOFF)
#################################################################
# Process negative control results into final table
#################################################################
n_false_rejections <- undercover_res |>
  filter(!(Method %in% c(c("NB regression (no covariates)", 
                           "NB regression (w/ covariates)", 
                           "SCEPTRE (no covariates)")))) |>
  group_by(dataset_rename, Method) |>
  summarize(n_false_reject = sum(p_value < alpha/n()),
                   Method = Method[1],
            `NT pairs` = n(), 
            .groups = "drop")

n_false_rejections_tab <- n_false_rejections |>
  pivot_wider(names_from = Method, values_from = n_false_reject) |>
  relocate("SCEPTRE", .after = "dataset_rename") |>
  relocate(`NT pairs`, .after = `KS test`) |>
  rename(Dataset = dataset_rename)

n_false_rejections_tab <- n_false_rejections_tab |>
  mutate(across(everything(), as.character)) |>
  bind_rows(
    n_false_rejections_tab |>
      summarise(across(-c(Dataset, `NT pairs`), mean)) |>
      mutate(Dataset = "Average") |>
      mutate(across(-Dataset, function(x)(as.character(round(x, 1)))))
  )  |>
  mutate(`NT pairs` = ifelse(is.na(`NT pairs`), "", `NT pairs`))
 n_false_rejections_tab <- n_false_rejections_tab[,c("Dataset", "SCEPTRE", "Seurat-Wilcox", "Seurat-NB", "t-test", "MAST", "KS test", "MIMOSCA", "NT pairs")]
 colnames(n_false_rejections_tab)[colnames(n_false_rejections_tab) == "Seurat-Wilcox"] <- "Seurat-\nWilcox"
 colnames(n_false_rejections_tab)[colnames(n_false_rejections_tab) == "Seurat-NB"] <- "Seurat-\nNB"
 colnames(n_false_rejections_tab)[colnames(n_false_rejections_tab) == "NT pairs"] <- "NT\npairs"
 
#################################################################
# Process positive control results
#################################################################

n_true_rejections_tab <- pc_res |>
  filter(n_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_control >= N_NONZERO_CONTROL_CUTOFF) |>
  group_by(dataset_rename, Method) |>
  summarize(n_pc_reject = sum(p_value < reject_thresh),
            `PC pairs` = n(),
            Method = Method[1], 
            .groups = "drop") |>
  group_by(dataset_rename) |>
  left_join(n_false_rejections,
            by = c("dataset_rename", "Method")) |>
  mutate(n_pc_reject = ifelse(n_false_reject <= max_false_reject, 
                              as.character(n_pc_reject), 
                              "-")) |>
  select(dataset_rename, Method, n_pc_reject, `PC pairs`) |>
  pivot_wider(names_from = Method, values_from = n_pc_reject) |>
  relocate("SCEPTRE", .after = "dataset_rename") |>
  relocate(`PC pairs`, .after = `KS test`) |>
  rename(Dataset = dataset_rename)
 
 n_true_rejections_tab <- n_true_rejections_tab[,c("Dataset", "SCEPTRE", "Seurat-Wilcox", "Seurat-NB", "t-test", "MAST", "KS test", "MIMOSCA", "PC pairs")]
 colnames(n_true_rejections_tab)[colnames(n_true_rejections_tab) == "Seurat-Wilcox"] <- "Seurat-\nWilcox"
 colnames(n_true_rejections_tab)[colnames(n_true_rejections_tab) == "Seurat-NB"] <- "Seurat-\nNB"
 colnames(n_true_rejections_tab)[colnames(n_true_rejections_tab) == "PC pairs"] <- "PC\npairs"
 
#################################################################
# Helper code for tables
#################################################################

# From https://cran.r-project.org/web/packages/gridExtra/vignettes/tableGrob.html
find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}

# set the theme for the table
table_theme <- ttheme_default(core=list(fg_params=list(hjust=1, x=0.9)),
                      base_size = 10)


#################################################################
# Format Type-I error table
#################################################################

# create gtable for negative control table
nt_table_g <- tableGrob(n_false_rejections_tab, theme = table_theme, rows = NULL)
  
# add bold font for lowest numbers of false rejections
for (row_idx in seq(1, nrow(n_false_rejections_tab))) {
  row <- as.integer(n_false_rejections_tab[row_idx, seq(2, 8L)])
  col_idxs <- which(row == min(row))
  for (col_idx in col_idxs) {
    nt_table_g$grobs[find_cell(nt_table_g, row_idx + 1L, col_idx + 1L, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
  }
}

# add horizontal line separating the average row from the rest
nt_table_g <- gtable_add_grob(nt_table_g,
                              grobs = segmentsGrob(
                                x0 = unit(0,"npc"),
                                y0 = unit(0,"npc"),
                                x1 = unit(1,"npc"),
                                y1 = unit(0,"npc"),
                                gp = gpar(lwd = 4.0)),
                              t = 8, b = 8, l = 1, r = 9)

# add vertical line separating the NT pairs column from the rest
nt_table_g <- gtable_add_grob(nt_table_g,
                              grobs = segmentsGrob( 
                                x0 = unit(0,"npc"),
                                y0 = unit(0,"npc"),
                                x1 = unit(0,"npc"),
                                y1 = unit(1,"npc"),
                                gp = gpar(lwd = 4.0)),
                              t = 8, b = 1, l = 9, r = 9)

# add vertical line separating the dataset column from the rest
nt_table_g <- gtable_add_grob(nt_table_g,
                              grobs = segmentsGrob( 
                                x0 = unit(0,"npc"),
                                y0 = unit(0,"npc"),
                                x1 = unit(0,"npc"),
                                y1 = unit(1,"npc"),
                                gp = gpar(lwd = 4.0)),
                              t = 9, b = 1, l = 2, r = 2)

# add title to the table
title <- textGrob("Number of false positives",gp=gpar(fontsize=12))
padding <- unit(5,"mm")
nt_table_g <- gtable_add_rows(
  nt_table_g, 
  heights = grobHeight(title) + padding,
  pos = 0)
nt_table_g <- gtable_add_grob(
  nt_table_g, 
  title, 
  1, 1, 1, ncol(nt_table_g))

# left-justify the first column
id <- which(grepl("core-fg", nt_table_g$layout$name ) & nt_table_g$layout$l == 1 )
for (i in id) {
  nt_table_g$grobs[[i]]$x <- unit(0.05, "npc")
  nt_table_g$grobs[[i]]$hjust <- 0
}

#################################################################
# Format power table
#################################################################

# create gtable for positive control table
pc_table_g <- tableGrob(n_true_rejections_tab, theme = table_theme, rows = NULL)

# add bold font for lowest numbers of false rejections
for (row_idx in seq(1, nrow(n_true_rejections_tab))) {
  row <- as.integer(n_true_rejections_tab[row_idx, seq(2, 8L)])
  row[is.na(row)] <- 0L
  col_idxs <- which(row == max(row))
  for (col_idx in col_idxs) {
    pc_table_g$grobs[find_cell(pc_table_g, row_idx + 1L, col_idx + 1L, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
  }
}

# add line separating the PC pairs column from the rest
pc_table_g <- gtable_add_grob(pc_table_g,
                grobs = segmentsGrob( 
                  x0 = unit(0,"npc"),
                  y0 = unit(0,"npc"),
                  x1 = unit(0,"npc"),
                  y1 = unit(1,"npc"),
                  gp = gpar(lwd = 4.0)),
                t = 7, b = 1, l = 9, r = 9)

# add line separating the dataset column from the rest
pc_table_g <- gtable_add_grob(pc_table_g,
                              grobs = segmentsGrob( 
                                x0 = unit(0,"npc"),
                                y0 = unit(0,"npc"),
                                x1 = unit(0,"npc"),
                                y1 = unit(1,"npc"),
                                gp = gpar(lwd = 4.0)),
                              t = 7, b = 1, l = 2, r = 2)

# add title
title <- textGrob("Number of true positives",gp=gpar(fontsize=12))
padding <- unit(5,"mm")
pc_table_g <- gtable_add_rows(
  pc_table_g, 
  heights = grobHeight(title) + padding,
  pos = 0)
pc_table_g <- gtable_add_grob(
  pc_table_g, 
  title, 
  1, 1, 1, ncol(pc_table_g))

# left justify the first column
id <- which(grepl("core-fg", pc_table_g$layout$name ) & pc_table_g$layout$l == 1 )
for (i in id) {
  pc_table_g$grobs[[i]]$x <- unit(0.05, "npc")
  pc_table_g$grobs[[i]]$hjust <- 0
}

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
          legend.position = c(0.25, 0.86),
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
final_plot <- ggarrange(
  ggarrange(qq_frangieh, qq_papalexi, nrow = 1),
  as_ggplot(nt_table_g),
  as_ggplot(pc_table_g),
  labels = "auto", 
  heights = c(1.2, 1, 0.8),
  ncol = 1
)

# define the file path
fig_4_filename <- paste0(
  .get_config_path("LOCAL_CODE_DIR"),
  "sceptre2-manuscript/R_scripts/figure_creation/fig_4/fig_4.png"
)

# save the figure
ggsave(filename = fig_4_filename, 
       plot = final_plot,
       device = "png", 
       width = 7,
       height = 8.5,
       bg = "white")
