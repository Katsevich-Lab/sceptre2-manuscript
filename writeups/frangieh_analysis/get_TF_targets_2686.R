###############################################################
# Description: Get the targets of each TF in the 2686 cell line 
# for the Frangieh enrichment analysis, based on intersecting 
# ATAC-seq data with transcription factor binding sites.
# 
# Author: Gene 
# 
# Date: April 21, 2023
###############################################################

###########################################################
# load libraries and resolve conflicts
###########################################################
library(GenomicRanges)
library(dplyr)
library(plyranges)
library(readr)
library(biomaRt)
library(conflicted)
library(data.table)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(base::intersect)

###########################################################
# set analysis parameters
###########################################################
JASPAR_THRESH <- 0.75
ATAC_THRESH <- 0.75
PROMOTER_WINDOW_WIDTH <- 5e3

###########################################################
# set filepaths
###########################################################
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
atac_filename <- "GSE205033_allpeaks_read.counts.rpkm.threshold.csv"
atac_data_fp <- paste0(sceptre2_dir, "/data/atacseq/", atac_filename)
method_results_fp <- paste0(
  sceptre2_dir, 
  "results/discovery_analyses/discovery_results_0423_processed.rds"
)
jaspar_tf_info_fp <- paste0(sceptre2_dir, "/data/jaspar/jaspar_tf_info.rds")
output_fp <- paste0(
  sceptre2_dir, 
  "results/discovery_analyses/TF_targets_2686.rds"
)

###########################################################
# read in the data and results
###########################################################
method_results <- readRDS(method_results_fp)
atac_data <- read_csv(atac_data_fp, show_col_types = FALSE, progress = FALSE)
jaspar_tf_info <- readRDS(jaspar_tf_info_fp)

###########################################################
# get TSS region for each gene
###########################################################

# get names of genes from Frangieh results
gene_names <- method_results |> 
  filter(dataset_rename == "Frangieh Control Gene") |> 
  pull(response_id) |>
  unique() |> 
  as.character()

# pull gene annotation info from ENSEMBL
ensembl <- useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', 
                      dataset = "hsapiens_gene_ensembl")
TSS_info <-getBM(attributes=c("external_gene_name", 
                              "chromosome_name", 
                              "start_position", 
                              "end_position", 
                              "strand"),
                 filters=c('external_gene_name'),
                 value = gene_names, mart=ensembl) |>
  filter(chromosome_name %in% c(1:22, "X", "Y")) |>
  mutate(TSS = ifelse(strand == 1, start_position, end_position),
         chromosome_name = paste0("chr", chromosome_name)) |>
  rename(gene = external_gene_name, chr = chromosome_name) |>
  select(gene, chr, TSS) |>
  group_by(gene, chr) |>
  summarise(TSS = round(mean(TSS)), .groups = "drop")

# save as GRanges object
TSS_GR <- GRanges(
  seqnames = TSS_info$chr,
  ranges = IRanges(start = TSS_info$TSS-PROMOTER_WINDOW_WIDTH, 
                   end = TSS_info$TSS+PROMOTER_WINDOW_WIDTH),
  gene = TSS_info$gene,
  TSS = TSS_info$TSS)

###########################################################
# get ATAC-seq peaks
###########################################################

# pull out the information for the 2686 cell type
atac_data_2686 <- atac_data |>
  # average biological replicates
  mutate(P2686 = (P2686A + P2686B + P2686C)/3) |> 
  filter(P2686 > quantile(P2686, ATAC_THRESH))

# save as GRanges object
atac_GR <- GRanges(
  seqnames = atac_data_2686$chrom,
  ranges = IRanges(start = atac_data_2686$chromStart, 
                   end = atac_data_2686$chromEnd)
)

###########################################################
# set up data structure for results
###########################################################

# get TFs for Frangieh data
TFs <- method_results |> 
  filter(dataset_rename == "Frangieh Control Gene") |> 
  pull(grna_group) |>
  intersect(jaspar_tf_info$name)

# get genes
genes <- TSS_GR$gene

# set up matrix to store results
target_info <- matrix(FALSE, 
                      nrow = length(genes), 
                      ncol = length(TFs),
                      dimnames = list(gene = genes, TF = TFs))

###########################################################
# get targets for each TF
###########################################################
for(TF in TFs){
  cat(sprintf("Working on %s...\n", TF))
  
  # get the motifs corresponding to this TF
  # (might have more than one motif per TF; use binding sites of all motifs)
  matrix_ids <- jaspar_tf_info |> 
    filter(name == TF) |>
    pull(matrix_id)
  
  # read the JASPAR locations of motifs
  jaspar_data <- lapply(matrix_ids, function(matrix_id){
    jaspar_filename <- paste0(matrix_id, ".tsv")
    jaspar_data_fp <- paste0(sceptre2_dir, "/data/jaspar/", jaspar_filename)
    jaspar_data <- read_tsv(jaspar_data_fp, 
                            col_names = c("chr", "start", "end", "TF", 
                                          "score1", "score2", "strand"), 
                            progress = FALSE,
                            show_col_types = FALSE)
  }) |> 
    rbindlist() |>
    filter(score1 > quantile(score1, JASPAR_THRESH))

  # convert to GRanges object  
  jaspar_GR <- GRanges(
    seqnames = jaspar_data$chr,
    ranges = IRanges(start = jaspar_data$start, end = jaspar_data$end)
  )
  
  # find ATAC-seq peaks overlapping a JASPAR binding site for this TF
  atac_overlapping_motif_GR <- atac_GR |>
    filter_by_overlaps(jaspar_GR)
  
  # find target genes as those where TSS region overlaps with TF binding site
  target_genes <- TSS_GR |> 
    filter_by_overlaps(atac_overlapping_motif_GR) |>
    as.data.frame() |>
    pull(gene)
  
  # set corresponding entries of output matrix to TRUE
  target_info[target_genes, TF] <- TRUE
}

# tidy and save
target_info |> 
  reshape2::melt(value.name = "target") |>
  as_tibble() |>
  saveRDS(output_fp)
cat(sprintf("Done!\n"))