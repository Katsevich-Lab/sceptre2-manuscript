# This script performs PCA on the datasets that exhibit miscalibration
#
# The Schraivogel TAP-seq and perturb-seq datasets exhibit some amount of miscalibration when
# analyzed using Fisher's exact test. As Fisher's exact test is an exact test of independence,
# this suggests that these datasets might be "contaminated" by latent confounding.

# This script performs PCA on the TAP-seq and perturb-seq datasets. It then (i)
# adds the top PC to the mattrix of cell-specific covariates and (ii) updates
# the formula object of the gene odm to include the top PC.
library(ondisc)
library(Seurat)

datasets <- c("ground_truth_tapseq", "ground_truth_perturbseq")
dataset_fps <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/schraivogel/", datasets)

########################
# 1. PCA on TAP-seq data
########################
dataset_fp <- dataset_fps[1]
gene_metadata_fp <- paste0(dataset_fp, "/gene/metadata_qc.rds")
gene_odm_fp <- paste0(dataset_fp, "/gene/matrix.odm")
odm <- read_odm(gene_odm_fp, gene_metadata_fp)
gene_mat <- as.matrix(lowmoi::load_whole_odm(odm))

# create seurat object
obj <- CreateSeuratObject(counts = as.matrix(lowmoi::load_whole_odm(odm)))
obj <- NormalizeData(obj)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = all.genes)
pc1_seurat <- obj@reductions$pca@cell.embeddings[,1]* -1
sds <- obj@reductions$pca@stdev
var_exp <- sds^2/sum(sds^2)
plot(var_exp) # the first PC explains >24% of the variance!
names(pc1_seurat) <- FALSE
hist(pc1_seurat)

# add leading pc as covariate
odm <- odm |> mutate_cell_covariates(pc_1 = pc1_seurat)
odm@misc$mimosca_formula <- update(odm@misc$mimosca_formula, ~ . + pc_1)
odm@misc$nb_regression_formula <- paste0(odm@misc$nb_regression_formula,  " + pc_1")
odm@misc$sceptre_formula <- update(odm@misc$sceptre_formula, ~ . + pc_1)
save_odm(odm = odm, metadata_fp = gene_metadata_fp)


#######################
# 2. PCA on perturb-seq
#######################
dataset_fp <- dataset_fps[2]
gene_metadata_fp <- paste0(dataset_fp, "/gene/metadata_qc.rds")
gene_odm_fp <- paste0(dataset_fp, "/gene/matrix.odm")
odm <- read_odm(gene_odm_fp, gene_metadata_fp)
gene_mat <- as.matrix(lowmoi::load_whole_odm(odm))

# create seurat object
obj <- CreateSeuratObject(counts = gene_mat)
obj <- NormalizeData(obj)
all.genes <- rownames(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
pc1_seurat <- obj@reductions$pca@cell.embeddings[,1]* -1
pc2_seurat <- obj@reductions$pca@cell.embeddings[,2]* -1
sds <- obj@reductions$pca@stdev
var_exp <- 100 * sds^2/sum(sds^2)
plot(pc1_seurat, pc2_seurat)
plot(var_exp) # the first two PCs explain > 28% of the variance
names(pc1_seurat) <- names(pc2_seurat) <- NULL

# Add leading PCs as covariates
odm <- odm |> mutate_cell_covariates(pc_1 = pc1_seurat, pc_2 = pc2_seurat)
odm@misc$mimosca_formula <- update(odm@misc$mimosca_formula, ~ . + pc_1 + pc_2)
odm@misc$nb_regression_formula <- paste0(odm@misc$nb_regression_formula,  " + pc_1 + pc_2")
odm@misc$sceptre_formula <- update(odm@misc$sceptre_formula, ~ . + pc_1 + pc_2)
save_odm(odm = odm, metadata_fp = gene_metadata_fp)
