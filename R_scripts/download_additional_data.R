library(tibble)
library(tidyr)
library(tidyverse)

sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

metadata_table <- 
  tribble(
    ~geo_id, ~TF, ~cell_type, ~condition, ~genome,~hTFtarget_id,~txt,
    935488,   "Stat1", "K562", "Ifng6h", "hg19",3385,F,
    935549,   "Irf1", "K562","Ifng6h","hg19",1762,F,
    1057025,"Irf1","Monocytes","Ifngh24h",'hg19',1767,T
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


