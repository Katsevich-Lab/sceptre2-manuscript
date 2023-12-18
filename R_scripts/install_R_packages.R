# print message
print("Installing required R packages...")

# pacman package for easy package handling
if (!require("pacman")) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
library(pacman)

# CRAN/Bioc
packs <- c("tidyverse", "readxl", "Seurat", "R.utils", "RCurl", "GenomicRanges",
           "MASS", "data.table", "rbioapi", "VGAM", "reticulate", "statmod",
           "SummarizedExperiment", "MAST")

for (pack in packs) p_install(pack)

# Github
github_packs <- c("katsevich-lab/lowmoi",
                  "timothy-barry/ondisc@5adcc53", # install commit 5adcc53 
                  "katsevich-lab/sceptre",
                  "katsevich-lab/katlabutils",
                  "katsevich-lab/sceptre",
                  "satijalab/seurat-data")

for (pack in github_packs) p_install_gh(pack)
