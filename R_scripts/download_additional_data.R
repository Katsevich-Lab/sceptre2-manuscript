library(tibble)
library(tidyr)
library(tidyverse)
library(readxl)
library(rbioapi)

sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

metadata_table <- 
  tribble(
    ~geo_id, ~TF, ~cell_type, ~condition, ~genome,~hTFtarget_id,~txt,
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
  filename <- paste0("GSM", geo_id, "_", toupper(TF), "peak_B.txt.gz")
  url <- paste0(geo_url, "GSM", geo_id, "&format=file&file=", filename)
  destfile <- paste0(chipseq_dir, "/", filename)
  download.file(url = url, destfile = destfile)
  R.utils::gunzip(destfile,overwrit = T)
}
