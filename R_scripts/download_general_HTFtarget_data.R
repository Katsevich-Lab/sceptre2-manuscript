
library(tibble)
library(tidyr)
library(tidyverse)

sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

metadata_table <- 
  tribble(
    ~TF,
    "STAT1",
    "IRF1",
    "STAT2",
    "SMAD4",
    "BRD4",
    "MYC"
  )


# directory for hTFtarget data
htftarget_dir <- paste0(sceptre2_dir, "/data/htftarget")
dir.create(htftarget_dir)



for(j in c(1:nrow(metadata_table))){
  TF <- metadata_table$TF[j]
  ############### Download hTFtarget data ###########################
  
  htftarget_url <- "http://bioinfo.life.hust.edu.cn/hTFtarget/static/hTFtarget/tmp_files/targets/"
  filename <- paste0(toupper(TF), ".target.txt.gz")
  destfile <- paste0(htftarget_dir, "/", filename)
  
  url <- paste0(htftarget_url, filename)
  download.file(url = url, destfile = destfile)
  R.utils::gunzip(destfile,overwrite = T)
  
}
