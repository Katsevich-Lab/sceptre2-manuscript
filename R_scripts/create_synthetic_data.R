library(ondisc)
set.seed(4)

# Randomly generate data: First, randomly sample gene-specific mean and size parameters from gamma distribution.
# Then, for each gene, randomly sample from an NB distribution.
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

# define hyperparameters
N_GENES <- 10000
N_GRNAS <- 35
N_CELLS <- 20000

# generate cell names, gene names, and gRNA names
cell_barcodes <- paste0("cell_", seq(1, N_CELLS))
gene_ids <- paste0("gene_", seq(1, N_GENES))
gRNA_ids <- paste0("NTC_", seq(1, N_GRNAS))

# create gene expression matrix
mus <- rgamma(n = N_GENES, shape = 1, rate = 2)
thetas <- runif(n = N_GENES, min = 5, max = 30)
gene_expression_mat <- sapply(X = seq(1, N_GENES), FUN = function(i) {
  MASS::rnegbin(N_CELLS, mus[i], thetas[i])
}) |> t()
rownames(gene_expression_mat) <- gene_ids
colnames(gene_expression_mat) <- cell_barcodes

# create gRNA expression matrix
gRNA_assignments <- sample(x = seq(1, N_GRNAS), size = N_CELLS, replace = TRUE)
gRNA_expression_mat <- sapply(X = gRNA_assignments, FUN = function(i) {
 out <- numeric(length = N_GRNAS)
 out[i] <- max(1, rpois(1, 100))
 return(out)
})
all(colSums(gRNA_expression_mat >= 1) == 1)
rownames(gRNA_expression_mat) <- gRNA_ids
colnames(gRNA_expression_mat) <- cell_barcodes

# perform quality control on the gene expression matrix
frac_cells_expressed <- rowMeans(gene_expression_mat >= 1)
gene_expression_mat <- gene_expression_mat[frac_cells_expressed > 0.005,]

# initialize ODM objects
to_save_fp_gene <- paste0(sceptre2_dir, "data/simulated/experiment_1/gene/matrix.odm")
create_ondisc_matrix_from_R_matrix(r_matrix = gene_expression_mat,
                                   barcodes = colnames(gene_expression_mat),
                                   features_df = data.frame(row.names(gene_expression_mat)),
                                   odm_fp = to_save_fp_gene,
                                   metadata_fp = paste0(sceptre2_dir, "data/simulated/experiment_1/gene/metadata_qc.rds"))

to_save_fp_gRNA <- paste0(sceptre2_dir, "data/simulated/experiment_1/grna/matrix.odm")
gRNA_odm <- create_ondisc_matrix_from_R_matrix(r_matrix = gRNA_expression_mat,
                                               barcodes = colnames(gRNA_expression_mat),
                                               features_df = data.frame(rownames(gRNA_expression_mat)),
                                               odm_fp = to_save_fp_gRNA)

# append target and target type to the gRNA odm
gRNA_odm <- gRNA_odm |>
  mutate_feature_covariates(target_type = "non-targeting", target = "non-targeting")
save_odm(odm = gRNA_odm, metadata_fp = paste0(sceptre2_dir, "data/simulated/experiment_1/grna/metadata_qc.rds"))
