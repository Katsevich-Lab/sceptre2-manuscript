# perform PCA on some subset of the datasets; add the top PCs to the matrix of technical factors
library(ondisc)

sceptre2_data_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")
tap_dir <- paste0(sceptre2_data_dir, "schraivogel/ground_truth_tapseq")
multimodal_fp <- paste0(tap_dir, "/multimodal_metadata.rds")
gene_fp <- paste0(tap_dir, "/gene/matrix.odm")
grna_assignment_fp <- paste0(tap_dir, "/grna_assignment/matrix.odm")

mm_odm <- read_multimodal_odm(odm_fps = c(gene_fp, grna_assignment_fp),
                    multimodal_metadata_fp = multimodal_fp)
# load gene expressions and normalize
gene_mat <- t(as.matrix(lowmoi::load_whole_odm(mm_odm@modalities$gene)))
lib_sizes <- rowSums(gene_mat)
gene_mat_norm <- log((gene_mat/lib_sizes) + 1)
# perform PCA
pca_fit <- prcomp(x = gene_mat_norm, center = TRUE, scale. = TRUE)
z <- pca_fit$x[,1]
names(z) <- NULL

# now, determine if pc_1 is related to gRNA presence/absence
grna_assignment_mat <- t(as.matrix(lowmoi::load_whole_odm(mm_odm@modalities$grna_assignment)))
p_vals <- numeric(length = ncol(grna_assignment_mat))
for (j in seq(1, ncol(grna_assignment_mat))) {
  print(j)
  x <- as.integer(grna_assignment_mat[,j])
  fit <- glm(formula = x ~ z, family = "binomial")
  s <- summary(fit)
  p_vals[j] <- s$coefficients[2, "Pr(>|z|)"]
}

# It appears to be the case that the first PC has a fairly strong association with many (but not necessarily all of) of the gRNAs.
# The first PC thus might be a confounder.
