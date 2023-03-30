sceptre2_dir <- .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")

metadata_table <- 
  tribble(
    ~geo_id, ~TF, ~cell_type, ~condition, hTFtarget_id,
    935488,   "Stat1", "K562", "Ifng6h", 3385
    "b",   2,
    "c",   3
  )

# directory for ChIP-seq data
chipseq_dir <- paste0(sceptre2_dir, "/data/chipseq")
dir.create(chipseq_dir)

# directory for hTFtarget data
htftarget_dir <- paste0(sceptre2_dir, "/data/htftarget")
dir.create(htftarget_dir)

geo_id <- 935488
build <- "hg19"
cell_type <- "K562"
TF <- "Stat1" # make a for loop to repeatedly download files for all TFs of interest
condition <- "Ifng6h"
hTFtarget_id <- 3385
geo_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc="

############### Download ENCODE ChIP-seq data ###########################

filename <- paste0("GSM", geo_id, 
                   "_", build, "_wgEncodeSydhTfbs", 
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
