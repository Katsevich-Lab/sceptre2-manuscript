library(tibble)
library(tidyr)
library(tidyverse)
library(readxl)
library(rbioapi)

sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

metadata_table <- 
  tribble(
    ~geo_id, ~TF, ~cell_type, ~condition, ~genome,~hTFtarget_id,~txt,
    935488,   "Stat1", "K562", "Ifng6h", "hg19",3385,F,
    935549,   "Irf1", "K562","Ifng6h","hg19",1762,F,
    1057025,"Irf1","Monocytes","Ifng24h",'hg19',1767,T,
    1057011,"Stat1","Monocytes","IFNg24h","hg19",3390,T
    )

# directory for ChIP-seq data
chipseq_dir <- paste0(sceptre2_dir, "/data/chipseq")
dir.create(chipseq_dir)

# directory for hTFtarget data
htftarget_dir <- paste0(sceptre2_dir, "/data/htftarget")
dir.create(htftarget_dir)

for(j in c(1:nrow(metadata_table))){
  geo_id <- metadata_table$geo_id[j]
  condition <- metadata_table$condition[j]
  cell_type <- metadata_table$cell_type[j]
  genome = metadata_table$genome[j]
  TF <- metadata_table$TF[j]
  condition <- metadata_table$condition[j]
  hTFtarget_id <- metadata_table$hTFtarget_id[j]
  txt = metadata_table$txt[j]
  
  geo_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc="
  
  ############### Download ENCODE ChIP-seq data ###########################
  if(txt == T){
    filename <- paste0("GSM", geo_id,
                       "_",toupper(TF),"peak_B.txt.gz" )
  }else{
    filename <- paste0("GSM", geo_id, 
                       "_", genome, "_wgEncodeSydhTfbs", 
                       cell_type, TF, condition, 
                       "StdPk.narrowPeak.gz")
  }
  
  url <- paste0(geo_url, "GSM", geo_id, "&format=file&file=", filename)
  destfile <- paste0(chipseq_dir, "/", filename)
  download.file(url = url, destfile = destfile)
  R.utils::gunzip(destfile,overwrit = T)
  
  ############### Download hTFtarget data ###########################
  
  htftarget_url <- "http://bioinfo.life.hust.edu.cn/hTFtarget/static/hTFtarget/tmp_files/targets/"
  filename <- paste0("dataset_", hTFtarget_id, ".", toupper(TF), ".target.txt.gz")
  destfile <- paste0(htftarget_dir, "/", filename)
  
  url <- paste0(htftarget_url, filename)
  download.file(url = url, destfile = destfile)
  R.utils::gunzip(destfile,overwrite = T)
  
}

############### Download ATAC-seq peaks for 2686 melanoma ###########################
atacseq_dir <- paste0(sceptre2_dir, "/data/atacseq")
dir.create(atacseq_dir)
filename <- "GSE205033_allpeaks_read.counts.rpkm.threshold.csv.gz"
destfile <- paste0(atacseq_dir, "/", filename)
download.file(url = paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205033/suppl/", filename),
              destfile = destfile)
R.utils::gunzip(destfile, overwrite = T)


##### Download ATAC seq data (THP-1 Cells with LPS Stimulation)
filename <- "GSM4425563_ATAC-seq_THP1_PMA_ctrl_TLR4_1hr.bw"
destfile <- paste0(atacseq_dir, "/", filename)
fileurl <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4425563&format=file&file=",
  filename
)
download.file(url = fileurl, destfile = destfile)

############### Download JASPAR TF binding sites ###########################

# create directory for JASPAR TF data
jaspar_dir <- paste0(sceptre2_dir, "/data/jaspar")
dir.create(jaspar_dir)

# download mapping between TFs and matrix IDs
jaspar_tf_info <- lapply(
  1:3,
  function(page) (rba_jaspar_collections_matrices(collection = "CORE", 
                                                  page = as.numeric(page), 
                                                  only_last_version = FALSE) %>% 
                    `$`(results) |> 
                    as_tibble()) 
) |>
  data.table::rbindlist() |>
  as_tibble() |>
  group_by(base_id, name) |>
  filter(version == max(version)) |>
  ungroup() |>
  select(matrix_id, name) |>
  mutate(name = case_when(
    name == "STAT1::STAT2" ~ "STAT2",
    name == "NFKB1" ~ "NFKBIA",
    .default = name
  ),
  matrix_id = ifelse(matrix_id == "MA0080.5", "MA0080.6", matrix_id))

saveRDS(jaspar_tf_info, paste0(jaspar_dir, "/jaspar_tf_info.rds"))

# download binding sites for TFs targeted in Frangieh 
# (note that for 2686 cells, we use hg38 rather than hg19)

frangieh_dir <- .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR")

# get list of genes targeted by Frangieh
frangieh_targets <- read_excel(
  path = paste0(frangieh_dir, "raw/supp_tables/41588_2021_779_MOESM3_ESM.xlsx"),
  sheet = 1,
  skip = 2
) |>
  rowwise() |>
  na.omit() |>
  filter(!grepl("SITE", `Guide Name`)) |>
  mutate(target = stringr::str_split_1(`Guide Name`, "_")[1]) |>
  pull(target) |>
  unique()

# get matrix IDs for TFs targeted by Papalexi
matrix_ids <- jaspar_tf_info |>
  filter(name %in% frangieh_targets) |>
  pull(matrix_id)

# for each matrix ID, download the corresponding binding sites
for (matrix_id in matrix_ids) {
  filename <- paste0(matrix_id, ".tsv.gz")
  destfile <- paste0(jaspar_dir, "/", filename)
  download.file(
    url = paste0(
      "http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/hg38/", 
      filename
    ),
    destfile = destfile
  )
  R.utils::gunzip(destfile, overwrite = T)
}

# download binding sites for TFs targeted in Papalexi 
# (note that for THP1 cells, we use hg19 rather than hg38)

papalexi_dir <- .get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")

# get list of genes targeted by Papalexi
papalexi_targets <- read_excel(
  path = paste0(papalexi_dir, "raw/41588_2021_778_MOESM4_ESM.xlsx"),
  sheet = 1
) |>
  mutate(`Target gene name` = ifelse(`Target gene name` == "NFKB1A", 
                                     "NFKBIA", 
                                     `Target gene name`)) |>
  pull(`Target gene name`)

# get matrix IDs for TFs targeted by Papalexi
matrix_ids <- jaspar_tf_info |>
  filter(name %in% papalexi_targets) |>
  pull(matrix_id)

# for each matrix ID, download the corresponding binding sites
for(matrix_id in matrix_ids){
  filename <- paste0(matrix_id, ".tsv.gz")
  filename_hg19 <- paste0(matrix_id,"_hg19", ".tsv.gz")
  destfile <- paste0(jaspar_dir, "/", filename_hg19)
  # for some reason the IRF1 binding sites are available in year 2020 but not 2022
  year <- if(matrix_id == "MA0050.2") 2020 else 2022 
  download.file(
    url = paste0(
      "http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/",
      year,
      "/hg19/",
      filename
    ),
    destfile = destfile
  )
  R.utils::gunzip(destfile, overwrite = T)  
}