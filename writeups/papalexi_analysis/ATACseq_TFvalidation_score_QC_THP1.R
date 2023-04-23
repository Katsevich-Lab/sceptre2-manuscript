library(GenomicRanges)
library(dplyr)
library(plyranges)
library(readr)
library(plotgardener)
library(GenomicRanges)
library(genomation)
library(biomaRt)
library(kableExtra)
library(varhandle)

#score QC to use
alpha = 0.50
sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")
data_dir <- sceptre2_dir
atac_dir <- paste0(sceptre2_dir, "/data/ATACseq")

# read in sceptre and seurat results
sceptre_path <- paste0(
  data_dir, "results/discovery_analyses/",
  "papalexi_gene_discovery_res.rds"
)
seurat_path <- paste0(
  data_dir, "results/papalexi_analysis/",
  "seurat_all_perturbations_results.rds"
)
sceptre_results <- readRDS(sceptre_path)
seurat_results <- readRDS(seurat_path)


# TSS information
gene_names <- sceptre_results |> 
  pull(response_id) |>
  unique() |> 
  as.character()

ensembl <- useEnsembl(host = 'https://grch37.ensembl.org',
                      biomart = 'ENSEMBL_MART_ENSEMBL', 
                      dataset = "hsapiens_gene_ensembl")



TSS_info <-getBM(attributes=c("external_gene_name", "chromosome_name", "start_position", 
                              "end_position", "strand"),
                 filters=c('external_gene_name'),
                 value = gene_names, mart=ensembl) |>
  filter(chromosome_name %in% c(1:22, "X", "Y")) |>
  mutate(TSS = ifelse(strand == 1, start_position, end_position),
         chromosome_name = paste0("chr", chromosome_name)) |>
  dplyr::rename(gene = external_gene_name, chr = chromosome_name) |>
  dplyr::select(gene, chr, TSS)



window_width <- 5e3
TSS_GR <- GRanges(
  seqnames = TSS_info$chr,
  ranges = IRanges(start = TSS_info$TSS-window_width, 
                   end = TSS_info$TSS+window_width),
  gene = TSS_info$gene,
  TSS = TSS_info$TSS)

# ATAC-seq data
#get filepath
atac_fp = paste0(atac_dir,"/GSM4425563_ATAC-seq_THP1_PMA_ctrl_TLR4_1hr.bw")
atac_data = readBigwig(
  file = paste0(atac_fp),
  chrom = NULL,
  chromstart = 1,
  chromend = .Machine$integer.max,
  strand = "*",
  params = NULL) |>
  filter(score > quantile(score,alpha)) |>
  filter(seqnames %in% c(as.character(1:22), "X", "Y"))%>%
  dplyr::rename(chrom = seqnames,chromStart = start, chromEnd = end)%>%
  mutate(chrom = paste0("chr", chrom))


#atac_data <- read_csv(atac_data_fp) |>
  #filter(P2686A > quantile(P2686A, alpha))
atac_GR <- GRanges(
  seqnames = atac_data$chrom,
  ranges = IRanges(start = atac_data$chromStart, end = atac_data$chromEnd)
)

# TF motif data
#get gene list
genes = as.data.frame(TSS_GR)
genes = unique(genes$gene)
#read jaspar association table
jaspar_tf_info <- readRDS(paste0(sceptre2_dir, "/data/jaspar/jaspar_tf_info_papalexi.rds"))
TFs <- jaspar_tf_info |> pull(name) |> unique()
#initialize empty matrix
TF_targets = matrix(0,length(genes),length(TFs))
colnames(TF_targets) = TFs
rownames(TF_targets) = genes
TF_targets = as.data.frame(TF_targets)
#iterate over each TF
for(TF in TFs){
  matrix_ids <- jaspar_tf_info |> 
    filter(name == TF) |>
    pull(matrix_id)
  
  # NOTE: Might have more than one motif per TF; use binding sites of all motifs
  jaspar_data <- lapply(matrix_ids, function(matrix_id){
    jaspar_filename <- paste0(matrix_id,"_hg19", ".tsv")
    jaspar_data_fp <- paste0(sceptre2_dir, "/data/jaspar/", jaspar_filename)
    jaspar_data <- read_tsv(jaspar_data_fp, 
                            col_names = c("chr", "start", "end", "TF", 
                                          "score1", "score2", "strand"))
  }) |> 
    data.table::rbindlist() |>
    filter(score1 > quantile(score1, alpha))
  
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
  
  TF_targets[[TF]] = as.integer(rownames(TF_targets)%in%target_genes)
  cat(sprintf("%s targets %s genes.\n", TF, length(target_genes)))
}

#save targets as rds file
target_fp = paste0(
  data_dir, "results/papalexi_analysis/",
  "ATACseq_TF_targets_score_QC_",100*alpha,".rds")

saveRDS(TF_targets,target_fp)
