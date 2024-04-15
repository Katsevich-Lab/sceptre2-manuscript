###############################################################
# Description: Get the targets of STAT1 and IRF1 in monocytes
# treated with IFN-gamma for the Papalexi enrichment analysis, 
# based on ChIP-seq data.
# 
# Author: Gene 
# 
# Date: April 22, 2023
###############################################################

###########################################################
# load libraries and resolve conflicts
###########################################################
library(GenomicRanges)
library(dplyr)
library(reshape2)
library(plyranges)
library(readr)
library(biomaRt)
library(data.table)
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(reshape2::melt)

###########################################################
# set analysis parameters
###########################################################
CHIPSEQ_THRESH <- 0.75
PROMOTER_WINDOW_WIDTH <- 5e3

###########################################################
# set filepaths
###########################################################
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
stat1_chipseq_filename <- "GSM1057011_STAT1peak_B.txt"
irf1_chipseq_filename <- "GSM1057025_IRF1peak_B.txt"
stat1_chipseq_fp <- paste0(sceptre2_dir, "/data/chipseq/", stat1_chipseq_filename)
irf1_chipseq_fp <- paste0(sceptre2_dir, "/data/chipseq/", irf1_chipseq_filename)
method_results_fp <- paste0(
  sceptre2_dir, 
  "results/discovery_analyses/discovery_results_0423_processed.rds"
)
output_fp <- paste0(
  sceptre2_dir, 
  "results/discovery_analyses/TF_targets_papalexi_chipseq.rds"
)

#######################################################################
# read in and process ChIP-seq data
#######################################################################

# ChIP-seq column names
chipseq_colnames <- c("chrom",'start_pos','end_pos',"pval","score",
  "pos_max_peak","max_peak_height", "rel_pos_max_peak_height",
  "peak_size","mid_point", "peak_to_mid_dist")

# read STAT1 ChIP-seq data
stat1_chipseq_data <- read_table(
  stat1_chipseq_fp,
  col_names = chipseq_colnames,
  show_col_types = FALSE
) |>
  filter(score > quantile(score, CHIPSEQ_THRESH)) |>
  mutate(TF = "STAT1")

# read IRF1 ChIP-seq data
irf1_chipseq_data <- read_table(
  irf1_chipseq_fp,
  col_names = chipseq_colnames,
  show_col_types = FALSE
) |>
  suppressWarnings() |>
  # change NA p-values to 1
  mutate(pval = ifelse(is.na(pval), 1, pval)) |>
  filter(score > quantile(score, CHIPSEQ_THRESH)) |>
  mutate(TF = "IRF1")

# combine the two into one data frame for convenience
chipseq_data <- bind_rows(stat1_chipseq_data |> 
                            select(chrom, start_pos, end_pos, score, TF), 
                          irf1_chipseq_data |>
                            select(chrom, start_pos, end_pos, score, TF))

# convert to GRanges object
chipseq_GR <- GRanges(
  seqnames = chipseq_data$chrom,
  ranges = IRanges(start = chipseq_data$start_pos, 
                   end = chipseq_data$end_pos),
  score = chipseq_data$score,
  TF = chipseq_data$TF)

###########################################################
# get TSS region for each gene
###########################################################

# read in method results
method_results <- readRDS(method_results_fp)

# get names of genes from Frangieh results
gene_names <- method_results |> 
  filter(dataset_rename == "Papalexi (Gene)") |> 
  pull(response_id) |>
  unique() |> 
  as.character()

# pull gene annotation info from ENSEMBL
ensembl <- useEnsembl(host = 'https://grch37.ensembl.org',
                      biomart = 'ENSEMBL_MART_ENSEMBL', 
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
  group_by(gene) |>
  # remove genes with multiple associated chromosomes
  filter(length(unique(chr)) == 1) |> 
  # for genes with multiple TSSs on same chromosome, average TSS positions
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
# get targets for each TF
###########################################################

# set up matrix to store results
TFs <- chipseq_data |> pull(TF) |> unique()
target_info <- matrix(FALSE, 
                      nrow = length(TSS_GR$gene), 
                      ncol = length(TFs),
                      dimnames = list(gene =  TSS_GR$gene, TF = TFs))

# loop over TFs
for(TF in TFs){
  cat(sprintf("Working on %s...\n", TF))
  
  # find target genes as those where TSS region overlaps with TF binding site
  target_genes <- TSS_GR |> 
    filter_by_overlaps(chipseq_GR |> plyranges::filter(TF == !!TF)) |>
    as.data.frame() |>
    pull(gene)
  
  # set corresponding entries of output matrix to TRUE
  target_info[target_genes, TF] <- TRUE
}

# tidy and save
target_info |> 
  melt(value.name = "target") |>
  as_tibble() |>
  saveRDS(output_fp)
cat(sprintf("Done!\n"))