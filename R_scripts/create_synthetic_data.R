library(ondisc)
set.seed(4)

# Randomly generate data: First, randomly sample gene-specific mean and size parameters from gamma distribution.
# Then, for each gene, randomly sample from an NB distribution.
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

# define hyperparameters
N_GENES <- 10000
N_GRNAS <- 35
N_NTC_GRNAS <- 30
N_CELLS <- 20000

# generate cell names, gene names, and grna names
cell_barcodes <- paste0("cell-", seq(1, N_CELLS))
gene_ids <- paste0("gene-", seq(1, N_GENES))
grna_ids <- c(paste0("NTC-", seq(1, N_NTC_GRNAS)), paste0("GENE-TARGET-", seq(1, N_GRNAS -  N_NTC_GRNAS)))

# create gene expression matrix
mus <- rgamma(n = N_GENES, shape = 1, rate = 2)
thetas <- runif(n = N_GENES, min = 5, max = 30)
gene_expression_mat <- sapply(X = seq(1, N_GENES), FUN = function(i) {
  MASS::rnegbin(N_CELLS, mus[i], thetas[i])
}) |> t()
rownames(gene_expression_mat) <- gene_ids
colnames(gene_expression_mat) <- cell_barcodes

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

# perform quality control on the gene expression matrix
# frac_cells_expressed <- rowMeans(gene_expression_mat >= 1)
# gene_expression_mat <- gene_expression_mat[frac_cells_expressed > 0.005,]

# initialize ODM objects
to_save_fp_gene <- paste0(sceptre2_dir, "data/simulated/experiment_1/gene/matrix.odm")
create_ondisc_matrix_from_R_matrix(r_matrix = gene_expression_mat,
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
