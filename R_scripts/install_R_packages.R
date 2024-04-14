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
                  "timothy-barry/ondisc@5adcc53", # install commit 5adcc53 (i.e, version 1.1.0)
                  "katsevich-lab/katlabutils",
                  "satijalab/seurat-data")
for (pack in github_packs) p_install_gh(pack)

# finally, install sceptre v0.3.0 (which is used in the benchmarking analysis and in figure 5)
devtools::install_github("katsevich-lab/sceptre", ref = "v0.3.0")
