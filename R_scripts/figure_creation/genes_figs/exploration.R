# Load packages
library(tidyverse)
library(katlabutils)
library(cowplot)
library(janitor)
library(kable)
library(kableExtra)
library(ggpubr)
library(grid)
library(gridExtra)
library(tableGrob)

# Load scripts and results
shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"), "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
result_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "results/")
undercover_res <- readRDS(paste0(result_dir,
                                 "undercover_grna_analysis/undercover_result_grp_1_processed.rds")) |>
  filter(n_nonzero_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_nonzero_control >= N_NONZERO_CONTROL_CUTOFF) |>
  mutate(Method = forcats::fct_recode(Method,
                                      "SCEPTRE" = "Sceptre",
                                      "Liscovitch" = "Liscovitch Method",
                                      "Schraivogel" = "Schraivogel Method",
                                      "Weissman" = "Weissman Method")) |>
  mutate(dataset_rename = forcats::fct_recode(dataset_rename,
                                              "Frangieh (Co Culture)" = "Frangieh Co Culture Gene",
                                              "Frangieh (Control)" = "Frangieh Control Gene",
                                              "Frangieh (IFN-\u03B3)" = "Frangieh Ifn Gamma Gene",
                                              "Papalexi (Gene)" = "Papalexi Eccite Screen Gene",
                                              "Papalexi (Protein)" = "Papalexi Eccite Screen Protein",
                                              "Schraivogel" = "Schraivogel Enhancer Screen",
                                              "Simulated" = "Simulated Experiment 1 Gene")) 
  
pc_res <- readRDS(paste0(result_dir, "positive_control_analysis/pc_results_processed.rds")) |>
  mutate(Method = forcats::fct_recode(Method,
                                      "SCEPTRE" = "Sceptre",
                                      "Liscovitch" = "Liscovitch Method",
                                      "Schraivogel" = "Schraivogel Method",
                                      "Weissman" = "Weissman Method")) |>
  mutate(dataset_rename = forcats::fct_recode(dataset_rename,
                                              "Frangieh (Co Culture)" = "Frangieh Co Culture Gene",
                                              "Frangieh (Control)" = "Frangieh Control Gene",
                                              "Frangieh (IFN-\u03B3)" = "Frangieh Ifn Gamma Gene",
                                              "Papalexi (Gene)" = "Papalexi Eccite Screen Gene",
                                              "Papalexi (Protein)" = "Papalexi Eccite Screen Protein",
                                              "Schraivogel" = "Schraivogel Enhancer Screen")) 
reject_thresh <- 1e-5


alpha <- 0.1
n_false_rejections <- undercover_res |>
  filter(!(Method %in% c(c("Nb Regression No Covariates", 
                           "Nb Regression W Covariates", 
                           "Sceptre No Covariates")))) |>
  group_by(dataset_rename, Method) |>
  summarize(n_false_reject = sum(p_value < alpha/n()),
                   Method = Method[1],
            `NT pairs` = n()) |>
  ungroup()

#################################################################
# Create tables
#################################################################

find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}

tt3 <- ttheme_default(core=list(fg_params=list(hjust=1, x=0.9)),
                      base_size = 10)


### Type-I error ###

n_false_rejections_tab <- n_false_rejections |>
  pivot_wider(names_from = Method, values_from = n_false_reject) |>
  relocate("SCEPTRE", .after = "dataset_rename") |>
  relocate(`NT pairs`, .after = `Weissman`) |>
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

nt_table_g <- tableGrob(n_false_rejections_tab, theme = tt3, rows = NULL)
  
nt_table_g$grobs[find_cell(nt_table_g, 2, 7, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 3, 4, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 3, 7, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 4, 7, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 5, 2, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 5, 4, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 6, 2, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 6, 4, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 6, 7, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 7, 5, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 8, 2, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 8, 6, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 8, 7, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
nt_table_g$grobs[find_cell(nt_table_g, 9, 2, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")


nt_table_g <- gtable_add_grob(nt_table_g,
                              grobs = segmentsGrob(
                                x0 = unit(0,"npc"),
                                y0 = unit(0,"npc"),
                                x1 = unit(1,"npc"),
                                y1 = unit(0,"npc"),
                                gp = gpar(lwd = 4.0)),
                              t = 8, b = 8, l = 1, r = 8)
nt_table_g <- gtable_add_grob(nt_table_g,
                              grobs = segmentsGrob( 
                                x0 = unit(0,"npc"),
                                y0 = unit(0,"npc"),
                                x1 = unit(0,"npc"),
                                y1 = unit(1,"npc"),
                                gp = gpar(lwd = 4.0)),
                              t = 8, b = 1, l = 8, r = 8)
nt_table_g <- gtable_add_grob(nt_table_g,
                              grobs = segmentsGrob( 
                                x0 = unit(0,"npc"),
                                y0 = unit(0,"npc"),
                                x1 = unit(0,"npc"),
                                y1 = unit(1,"npc"),
                                gp = gpar(lwd = 4.0)),
                              t = 9, b = 1, l = 2, r = 2)


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


id <- which(grepl("core-fg", nt_table_g$layout$name ) & nt_table_g$layout$l == 1 )

# loop through grobs and change relevant parts
for (i in id) {
  nt_table_g$grobs[[i]]$x <- unit(0.05, "npc")
  nt_table_g$grobs[[i]]$hjust <- 0
}
plot(nt_table_g)

### Power ###

n_pc_reject_df <- pc_res |>
  filter(n_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_control >= N_NONZERO_CONTROL_CUTOFF) |>
  group_by(dataset_rename, Method) |>
  summarize(n_pc_reject = sum(p_value < reject_thresh),
            `PC pairs` = n(),
            Method = Method[1]) |>
  group_by(dataset_rename) |>
  left_join(n_false_rejections,
            by = c("dataset_rename", "Method")) |>
  mutate(n_pc_reject = ifelse(n_false_reject <= 50, 
                              as.character(n_pc_reject), 
                              "-")) |>
  select(dataset_rename, Method, n_pc_reject, `PC pairs`) |>
  ungroup()

type_I_err_table_filename <- paste0(.get_config_path("LOCAL_CODE_DIR"), 
                           "sceptre2-manuscript/R_scripts/figure_creation/genes_figs/type_I_err_table.png")

pc_table <- n_pc_reject_df |>
  pivot_wider(names_from = Method, values_from = n_pc_reject) |>
  relocate("SCEPTRE", .after = "dataset_rename") |>
  relocate(`PC pairs`, .after = `Weissman`) |>
  rename(Dataset = dataset_rename)

pc_table_g <- tableGrob(pc_table, theme = tt3, rows = NULL)
pc_table_g$grobs[find_cell(pc_table_g, 2, 2, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g$grobs[find_cell(pc_table_g, 3, 2, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g$grobs[find_cell(pc_table_g, 4, 2, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g$grobs[find_cell(pc_table_g, 5, 2, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g$grobs[find_cell(pc_table_g, 6, 2, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g$grobs[find_cell(pc_table_g, 6, 3, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g$grobs[find_cell(pc_table_g, 6, 4, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g$grobs[find_cell(pc_table_g, 6, 5, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g$grobs[find_cell(pc_table_g, 6, 6, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g$grobs[find_cell(pc_table_g, 6, 7, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g$grobs[find_cell(pc_table_g, 7, 3, "core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
pc_table_g <- gtable_add_grob(pc_table_g,
                grobs = segmentsGrob( # line across the bottom
                  x0 = unit(0,"npc"),
                  y0 = unit(0,"npc"),
                  x1 = unit(0,"npc"),
                  y1 = unit(1,"npc"),
                  gp = gpar(lwd = 4.0)),
                t = 7, b = 1, l = 8, r = 8)
pc_table_g <- gtable_add_grob(pc_table_g,
                              grobs = segmentsGrob( 
                                x0 = unit(0,"npc"),
                                y0 = unit(0,"npc"),
                                x1 = unit(0,"npc"),
                                y1 = unit(1,"npc"),
                                gp = gpar(lwd = 4.0)),
                              t = 7, b = 1, l = 2, r = 2)

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

id <- which(grepl("core-fg", pc_table_g$layout$name ) & pc_table_g$layout$l == 1 )

# loop through grobs and change relevant parts
for (i in id) {
  pc_table_g$grobs[[i]]$x <- unit(0.05, "npc")
  pc_table_g$grobs[[i]]$hjust <- 0
}
plot(pc_table_g)

###########################

my_values <- my_cols[names(my_cols) %in% c("Seurat De", "SCEPTRE")]

qq_frangieh <- undercover_res |>
  mutate(Method = fct_recode(Method,
                             "SCEPTRE" = "Sceptre")) |>
  mutate(Method = fct_relevel(Method, "SCEPTRE", after = Inf)) |>
  filter(dataset == "frangieh_ifn_gamma_gene",
         method %in% c("sceptre", "seurat_de")) |> 
  ggplot(mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.85) +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10)) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
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

qq_papalexi <- undercover_res |>
  mutate(Method = fct_recode(Method,
                             "SCEPTRE" = "Sceptre")) |>
  mutate(Method = fct_relevel(Method, "SCEPTRE", after = Inf)) |>
  filter(dataset == "papalexi_eccite_screen_gene",
         method %in% c("sceptre", "seurat_de")) |> 
  ggplot(mapping = aes(y = p_value, col = Method)) +
  stat_qq_points(ymin = 1e-8, size = 0.85) +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10)) +
  labs(x = "Expected null p-value", y = "Observed p-value") +
  geom_abline(col = "black") +
  ggtitle("Papalexi (gene) neg. controls") +
  scale_color_manual(values = my_values) + 
  my_theme_no_legend



###########################

final_plot <- ggarrange(
  ggarrange(qq_frangieh, qq_papalexi, nrow = 1),
  as_ggplot(nt_table_g),
  as_ggplot(pc_table_g),
  labels = "auto", 
  heights = c(1.2, 1, 0.8),
  ncol = 1
)

final_plot

ggarrange(as_ggplot(nt_table_g), 
          as_ggplot(pc_table_g),
          labels = "auto",
          ncol = 1)



ggsave(filename = "~/Desktop/final_plot.png", 
       plot = final_plot, 
       device = "png", 
       width = 7.5, 
       height = 8.5)

#################################################################
# Proportion of rejections
#################################################################

#prop_reject_df <- 
  pc_res |>
  filter(n_treatment >= N_NONZERO_TREATMENT_CUTOFF,
         n_control >= N_NONZERO_CONTROL_CUTOFF) |>
  group_by(dataset_rename, Method) |>
  summarize(prop_pc_reject = mean(p_value < reject_thresh),
            Method = Method[1]) |>
  group_by(dataset_rename) |>
  left_join(n_false_rejections,
            by = c("dataset_rename", "Method")) |>
  filter(n_false_reject < 50) |>
  select(dataset_rename, Method, prop_pc_reject) |>
  ungroup() |>
  pivot_wider(names_from = "Method", values_from = "prop_pc_reject")

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