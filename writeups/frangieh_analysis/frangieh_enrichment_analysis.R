library(GenomicRanges)
library(plyranges)

sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
atac_filename <- "GSE205033_allpeaks_read.counts.rpkm.threshold.csv"
atac_data_fp <- paste0(sceptre2_dir, "/data/atacseq/", atac_filename)
sceptre_results_fp <- paste0(
  sceptre2_dir, 
  "results/discovery_analyses/frangieh_discovery_res.rds"
)

sceptre_results <- readRDS(sceptre_results_fp)

# TSS information
gene_names <- sceptre_results |> 
  pull(response_id) |>
  unique() |> 
  as.character()

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
  select(gene, chr, TSS)

window_width <- 5e3
TSS_GR <- GRanges(
  seqnames = TSS_info$chr,
  ranges = IRanges(start = TSS_info$TSS-window_width, 
                   end = TSS_info$TSS+window_width),
  gene = TSS_info$gene,
  TSS = TSS_info$TSS)

# ATAC-seq data
atac_data <- read_csv(atac_data_fp) |>
  filter(P2686A > quantile(P2686A, 0.75))
atac_GR <- GRanges(
  seqnames = atac_data$chrom,
  ranges = IRanges(start = atac_data$chromStart, end = atac_data$chromEnd)
) |>

# TF motif data
jaspar_tf_info <- readRDS(paste0(sceptre2_dir, "/data/jaspar/jaspar_tf_info.rds"))
TFs <- jaspar_tf_info |> pull(name) |> unique()

for(TF in TFs){
  matrix_ids <- jaspar_tf_info |> 
    filter(name == TF) |>
    pull(matrix_id)
  
  # NOTE: Might have more than one motif per TF; use binding sites of all motifs
  jaspar_data <- lapply(matrix_ids, function(matrix_id){
    jaspar_filename <- paste0(matrix_id, ".tsv")
    jaspar_data_fp <- paste0(sceptre2_dir, "/data/jaspar/", jaspar_filename)
    jaspar_data <- read_tsv(jaspar_data_fp, 
                            col_names = c("chr", "start", "end", "TF", 
                                          "score1", "score2", "strand"))
  }) |> 
    data.table::rbindlist() |>
    filter(score1 > quantile(score1, 0.75))
  
  jaspar_GR <- GRanges(
    seqnames = jaspar_data$chr,
    ranges = IRanges(start = jaspar_data$start, end = jaspar_data$end)
  )
  
  atac_overlapping_motif_GR <- atac_GR |>
    filter_by_overlaps(jaspar_GR)
  
  target_genes <- TSS_GR |> 
    filter_by_overlaps(atac_overlapping_motif_GR) |>
    as.data.frame() |>
    pull(gene)
  
  cat(sprintf("%s targets %s genes.\n", TF, length(target_genes)))
}