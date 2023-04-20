library(tibble)
library(tidyr)
library(tidyverse)

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

# directory for hTFtarget data
chromhmm_dir <- paste0(sceptre2_dir, "data/ChromHMM")
dir.create(chromhmm_dir)


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

############### Download ChromHMM data for Monocytes ###########################

chromhmm_url <- "https://www.encodeproject.org/files/"
ID = "ENCFF269WBG"
filename = paste0(ID,"/@@download/",ID,".bed.gz")
destfile <- paste0(chromhmm_dir, "/", ID,".bed.gz")
url <- paste0(chromhmm_url, filename)
download.file(url = url, destfile = destfile)
R.utils::gunzip(destfile,overwrite = T)

############### Download ChromHMM data for K562 ###########################

chromhmm_url <- "https://www.encodeproject.org/files/"
ID = "ENCFF163QUM"
filename = paste0(ID,"/@@download/",ID,".bed.gz")
destfile <- paste0(chromhmm_dir, "/", ID,".bed.gz")
url <- paste0(chromhmm_url, filename)
download.file(url = url, destfile = destfile)
R.utils::gunzip(destfile,overwrite = T)

# http://bioinfo.life.hust.edu.cn/hTFtarget/static/hTFtarget/tmp_files/targets/dataset_3390.STAT1.target.txt.gz
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1057011&format=file&file=GSM1057011%5FSTAT1peak%5FB%2Etxt%2Egz

############### Download ATAC-seq peaks for 2686 melanoma ###########################

# directory for ChIP-seq data
atacseq_dir <- paste0(sceptre2_dir, "/data/atacseq")
dir.create(atacseq_dir)
filename <- "GSE205033_allpeaks_read.counts.rpkm.threshold.csv.gz"
destfile <- paste0(atacseq_dir, "/", filename)
download.file(url = paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205033/suppl/", filename),
              destfile = destfile)
R.utils::gunzip(destfile, overwrite = T)

############### Download JASPAR TF binding sites ###########################

jaspar_dir <- paste0(sceptre2_dir, "/data/jaspar")
dir.create(jaspar_dir)

jaspar_tf_info <- lapply(
  1:2,
  function(page) (rba_jaspar_collections_matrices(collection = "CORE", 
                                                  page = as.numeric(page), 
                                                  only_last_version = TRUE) %>% 
                    `$`(results) |> 
                    as_tibble() |> 
                    select(matrix_id, name))
) |>
  data.table::rbindlist() |>
  as_tibble()

frangieh_dir <- .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR")
frangieh_targets <- readxl::read_excel(path = paste0(frangieh_dir, "raw/supp_tables/41588_2021_779_MOESM3_ESM.xlsx"), 
                                       sheet = 1,
                                       skip = 2) |>
  rowwise() |>
  na.omit() |>
  filter(!grepl("SITE", `Guide Name`)) |>
  mutate(target = stringr::str_split_1(`Guide Name`, "_")[1]) |>
  pull(target) |>
  unique()

matrix_ids <- jaspar_tf_info |>
  filter(name %in% frangieh_targets) |>
  pull(matrix_id)

jaspar_tf_info |>
  filter(name %in% frangieh_targets) |>
  saveRDS(paste0(jaspar_dir, "/jaspar_tf_info.rds"))
  
for(matrix_id in matrix_ids){
  filename <- paste0(matrix_id, ".tsv.gz")
  destfile <- paste0(jaspar_dir, "/", filename)
  download.file(url = paste0("http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/hg38/", filename),
                destfile = destfile)
  R.utils::gunzip(destfile, overwrite = T)  
}