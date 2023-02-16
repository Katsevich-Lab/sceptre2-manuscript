library(ondisc)
library(sceptre2)
library(sceptre)

#devtools::install_github('timothy-barry/ondisc')

LOCAL_SCEPTRE2_DATA_DIR <-.get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
papalexi_dir <- paste0(LOCAL_SCEPTRE2_DATA_DIR, "data/papalexi/eccite_screen/")

# gene info
gene_odm_fp <- paste0(papalexi_dir, "gene/matrix.odm")

# grna info
grna_odm_fp <- paste0(papalexi_dir, "grna_assignment/matrix.odm")

# protein info
protein_odm_fp <- paste0(papalexi_dir, "protein/matrix.odm")

# mm odm metadata fp
mm_metadata_fp <- paste0(papalexi_dir, "multimodal_metadata.rds")

# construct mm odm
mm_odm <- read_multimodal_odm(odm_fps = c(gene_odm_fp, grna_odm_fp, protein_odm_fp),
                              multimodal_metadata_fp = mm_metadata_fp)

# set formulas, grna group target name
gene_formula <- ~ log(gene_n_umis) + log(gene_n_nonzero) + bio_rep + phase + p_mito
protein_formula <- ~ log(protein_n_umis) + bio_rep + phase + p_mito
grna_group <- "target"

# set hyperparameters
B <- 2500000
side <- "both"
max_b_per_batch <- 250000
in_memory <- TRUE
statistic <- "full" # "distilled" is faster but might be less powerful
return_dist <- FALSE
screen_b <- 25000

# randomly select gene-grna group pairs to analyze
gene_grna_group_pairs <- expand.grid(response_id = mm_odm |>
                                       get_modality("gene") |>
                                       get_feature_ids(),
                                     grna_group = c("CUL3"))

grna_all = unique(mm_odm@modalities$grna_assignment@feature_covariates$target)[-c(3,4)]
# select protein-grna group pairs to analyze
#CD274 perturbed cells not present in sample
# table(substr(mm_odm@modalities$grna_assignment@cell_covariates$assigned_grna,1,4))
protein_grna_group_pairs <- expand.grid(response_id = mm_odm |>
                                          get_modality("protein") |>
                                          get_feature_ids(),
                                        grna_group = c("IRF1","BRD4","CUL3","CMTM6","ATF2",
                                                       "CAV1","CD86","ETV7","IFNGR1","IFNGR2","IRF7",
                                                       "JAK2","MARCH8","MYC","NFKBIA","PDCD1LG2",
                                                       "POU2F2","SMAD4","SPI1","STAT1","STAT2","STAT3",
                                                       "STAT5A","TNFRSF14","UBE2L6"))

# analyze the gene data
gene_result <- run_sceptre_low_moi(mm_odm = mm_odm,
                                   response_grna_group_pairs = gene_grna_group_pairs,
                                   form = gene_formula,
                                   response_modality_name = "gene",
                                   grna_modality_name = "grna_assignment",
                                   grna_group_column_name = "target",
                                   B = B,
                                   side = side,
                                   max_b_per_batch = max_b_per_batch,
                                   in_memory = in_memory,
                                   statistic = statistic,
                                   return_dist = return_dist,
                                   screen_b = screen_b)

# analyze the protein data
protein_result <- run_sceptre_low_moi(mm_odm,
                                      response_grna_group_pairs = protein_grna_group_pairs,
                                      form = protein_formula,
                                      response_modality_name = "protein",
                                      grna_modality_name = "grna_assignment",
                                      grna_group_column_name = "target",
                                      B = B,
                                      side = side,
                                      max_b_per_batch = max_b_per_batch,
                                      in_memory = in_memory,
                                      statistic  = statistic,
                                      return_dist = return_dist,
                                      screen_b = screen_b)

#save as RDS files. this code moves it to the wrong folders though
saveRDS(gene_result,'gene_result.rds')
saveRDS(protein_result,'protein_result.rds')


#get pvalues from sceptre
P_adj = protein_result[,1]
#unlist pvalues
P_adj = unlist(P_adj)
#some pvalues are negative so take absolute value
P_adj = abs(P_adj)
#make numeric 
P_adj = as.numeric(P_adj)
#perform BH procedure
P_adj = p.adjust(P_adj,method = 'BH')

#replace results matrix pvalues wioth adjusted pvalues
protein_adjusted= cbind(P_adj,protein_result[,c(2,3)])
#filter to just look at PDL1 pvalues
protein_adjusted = protein_adjusted[which(protein_adjusted[,3] == 'PDL1')]
#view
View(protein_adjusted)

#get pvalues from sceptre
P_adj = gene_result[,1]
#unlist pvalues
P_adj = unlist(P_adj)
#some pvalues are negative so take absolute value
P_adj = abs(P_adj)
#make numeric 
P_adj = as.numeric(P_adj)
#perform BH procedure
P_adj = p.adjust(P_adj,method = 'BH')

#get pvalue corresponding to PDL1
P_adj[which(gene_result[,3] == 'CD274')]

