library(ondisc)
set.seed(4)

###########
# DATASET 1
###########
# Randomly generate data: First, randomly sample gene-specific mean and size parameters from gamma distribution.
# Then, for each gene, randomly sample from an NB distribution.
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
# define hyperparameters
N_GENES <- 5000
N_GRNAS <- 30
N_NTC_GRNAS <- 25
N_CELLS <- 10000

# generate cell names, gene names, and grna names
cell_barcodes <- paste0("cell-", seq(1, N_CELLS))
gene_ids <- paste0("gene-", seq(1, N_GENES))
grna_ids <- c(paste0("NTC-", seq(1, N_NTC_GRNAS)),
              paste0("GENE-TARGET-", seq(1, N_GRNAS -  N_NTC_GRNAS)))

########################
# NEGATIVE CONTROL PAIRS
########################
# create grna expression matrix
grna_assignments <- sample(x = seq(1, N_GRNAS), size = N_CELLS, replace = TRUE)
grna_expression_mat <- sapply(X = grna_assignments, FUN = function(i) {
  out <- numeric(length = N_GRNAS)
  out[i] <- max(1, rpois(1, 100))
  return(out)
})
all(colSums(grna_expression_mat >= 1) == 1)
rownames(grna_expression_mat) <- grna_ids
colnames(grna_expression_mat) <- cell_barcodes

# create gene expression matrix
mus <- rgamma(n = N_GENES, shape = 0.5, rate = 2)
thetas <- runif(n = N_GENES, min = 1, max = 25)
gene_expression_mat <- sapply(X = seq(1, N_GENES), FUN = function(i) {
  MASS::rnegbin(N_CELLS, mus[i], thetas[i])
}) |> t()
rownames(gene_expression_mat) <- gene_ids
colnames(gene_expression_mat) <- cell_barcodes

# initialize ODM objects
to_save_fp_gene <- paste0(sceptre2_dir, "data/simulated/experiment_1/gene/matrix.odm")
gene_odm <- create_ondisc_matrix_from_R_matrix(r_matrix = gene_expression_mat,
                                               barcodes = colnames(gene_expression_mat),
                                               features_df = data.frame(row.names(gene_expression_mat)),
                                               odm_fp = to_save_fp_gene,
                                               metadata_fp = paste0(sceptre2_dir, "data/simulated/experiment_1/gene/metadata_orig.rds"))

to_save_fp_grna <- paste0(sceptre2_dir, "data/simulated/experiment_1/grna_expression/matrix.odm")
grna_odm <- create_ondisc_matrix_from_R_matrix(r_matrix = grna_expression_mat,
                                               barcodes = colnames(grna_expression_mat),
                                               features_df = data.frame(rownames(grna_expression_mat)),
                                               odm_fp = to_save_fp_grna)

# append target and target type to the grna odm
grna_tbl <- data.frame(target_type = c(rep("non-targeting", N_NTC_GRNAS), rep("gene", N_GRNAS - N_NTC_GRNAS)),
                       target = c(rep("non-targeting", N_NTC_GRNAS), rep("gene-1", N_GRNAS - N_NTC_GRNAS)))

# update the grna odm
grna_odm <- grna_odm |>
  mutate_feature_covariates(grna_tbl)
save_odm(odm = grna_odm, metadata_fp = paste0(sceptre2_dir, "data/simulated/experiment_1/grna_expression/metadata_orig.rds"))

# finally, create the matrix of grna assignments
convert_assign_list_to_sparse_odm(cell_barcodes = cell_barcodes,
                                  grna_ids = grna_ids,
                                  grna_assignment_list = as.list(grna_ids[grna_assignments]),
                                  odm_fp = paste0(sceptre2_dir, "data/simulated/experiment_1/grna_assignment/matrix.odm"),
                                  metadata_fp = paste0(sceptre2_dir, "data/simulated/experiment_1/grna_assignment/metadata_orig.rds"),
                                  features_metadata_df = grna_tbl)

########################
# POSITIVE CONTROL PAIRS
########################
# set the number of grnas
set.seed(4)
N_PC_GRNAS <- 25
N_NT_GRNAS <- 5
N_GRNAS <- N_PC_GRNAS + N_NT_GRNAS 
N_GENES <- N_GRNAS
N_CELLS <- 10000
grna_assignments <- sample(x = seq(1, N_GRNAS), size = N_CELLS, replace = TRUE)
mu_0s <- rgamma(n = N_GENES, shape = 2, rate = 0.5)
mu_1s <- rgamma(n = N_GENES, shape = 2, rate = 0.5)
thetas <- runif(n = N_GENES, min = 1, max = 25)

# generate the gene expression matrix
gene_expression_mat <- sapply(X = seq(1, N_GENES), FUN = function(i) {
  alt_cell_idxs <- which(grna_assignments == i)
  null_cell_idxs <- which(grna_assignments != i)
  y_null <- MASS::rnegbin(length(null_cell_idxs), mu_0s[i], thetas[i])
  y_alt <- MASS::rnegbin(length(alt_cell_idxs), mu_1s[i], thetas[i])
  y <- integer(N_CELLS)
  y[null_cell_idxs] <- y_null
  y[alt_cell_idxs] <- y_alt
  return(y)
}) |> t()

# generate the gRNA expression matrix
grna_expression_mat <- sapply(X = grna_assignments, FUN = function(i) {
  out <- numeric(length = N_GRNAS)
  out[i] <- max(1, rpois(1, 100))
  return(out)
})

# generate the gene ids, grna ids, and cell barcodes; append to the ondisc objects
cell_barcodes <- paste0("cell-", seq(1, N_CELLS))
gene_ids <- paste0("gene-", seq(1, N_GENES))
grna_ids <- c(paste0("GENE-TARGET-", seq(1, N_PC_GRNAS)),
              paste0("NTC-", seq(1, N_GRNAS - N_PC_GRNAS)))

rownames(gene_expression_mat) <- gene_ids
colnames(gene_expression_mat) <- cell_barcodes
rownames(grna_expression_mat) <- grna_ids
colnames(grna_expression_mat) <- cell_barcodes

# initialize ODM objects
# gene
to_save_fp_gene <- paste0(sceptre2_dir, "data/simulated/experiment_2/gene/matrix.odm")
gene_odm <- create_ondisc_matrix_from_R_matrix(r_matrix = gene_expression_mat,
                                               barcodes = colnames(gene_expression_mat),
                                               features_df = data.frame(row.names(gene_expression_mat)),
                                               odm_fp = to_save_fp_gene,
                                               metadata_fp = paste0(sceptre2_dir, "data/simulated/experiment_2/gene/metadata_orig.rds"))

# grna
to_save_fp_grna <- paste0(sceptre2_dir, "data/simulated/experiment_2/grna_expression/matrix.odm")
grna_odm <- create_ondisc_matrix_from_R_matrix(r_matrix = grna_expression_mat,
                                               barcodes = colnames(grna_expression_mat),
                                               features_df = data.frame(rownames(grna_expression_mat)),
                                               odm_fp = to_save_fp_grna)

# append target and target type to the grna odm
grna_tbl <- data.frame(target_type = c(rep("gene", N_PC_GRNAS), rep("non-targeting", N_NT_GRNAS)),
                       target = c(paste0("gene-", seq(1, N_PC_GRNAS)), rep("non-targeting", N_NT_GRNAS)))

# update the grna odm
grna_odm <- grna_odm |>
  mutate_feature_covariates(grna_tbl)
save_odm(odm = grna_odm, metadata_fp = paste0(sceptre2_dir, "data/simulated/experiment_2/grna_expression/metadata_orig.rds"))

# finally, create the matrix of grna assignments
convert_assign_list_to_sparse_odm(cell_barcodes = cell_barcodes,
                                  grna_ids = grna_ids,
                                  grna_assignment_list = as.list(grna_ids[grna_assignments]),
                                  odm_fp = paste0(sceptre2_dir, "data/simulated/experiment_2/grna_assignment/matrix.odm"),
                                  metadata_fp = paste0(sceptre2_dir, "data/simulated/experiment_2/grna_assignment/metadata_orig.rds"),
                                  features_metadata_df = grna_tbl)
