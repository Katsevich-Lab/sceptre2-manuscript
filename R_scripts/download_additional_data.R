library(tibble)
library(tidyr)
library(tidyverse)

sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

metadata_table <- 
  tribble(
    ~geo_id, ~TF, ~cell_type, ~condition, ~genome,~hTFtarget_id,
    935488,   "Stat1", "K562", "Ifng6h", "hg19",3385,
    935549,   "Irf1", "K562","Ifng6h","hg19",1762
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
  
  geo_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc="
  
  ############### Download ENCODE ChIP-seq data ###########################
  
  filename <- paste0("GSM", geo_id, 
                     "_", genome, "_wgEncodeSydhTfbs", 
                     cell_type, TF, condition, 
                     "StdPk.narrowPeak.gz")
  url <- paste0(geo_url, "GSM", geo_id, "&format=file&file=", filename)
  destfile <- paste0(chipseq_dir, "/", filename)
  download.file(url = url, destfile = destfile)
  R.utils::gunzip(destfile)
  
  ############### Download hTFtarget data ###########################
  
  htftarget_url <- "http://bioinfo.life.hust.edu.cn/hTFtarget/static/hTFtarget/tmp_files/targets/"
  filename <- paste0("dataset_", hTFtarget_id, ".", toupper(TF), ".target.txt.gz")
  destfile <- paste0(htftarget_dir, "/", filename)
  
  url <- paste0(htftarget_url, filename)
  download.file(url = url, destfile = destfile)
  R.utils::gunzip(destfile)
  
}




