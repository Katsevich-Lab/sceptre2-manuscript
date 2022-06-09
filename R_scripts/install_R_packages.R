# pacman package for easy package handling
if (!require("pacman")) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
library(pacman)

# CRAN
p_install(readr, force = FALSE)
p_install(dplyr, force = FALSE)
p_install(readxl, force = FALSE)
p_install(Seurat, force = FALSE)
p_install(R.utils, force = FALSE)
p_install(purrr, force = FALSE)
p_install(RCurl, force = FALSE)
p_install(GenomicRanges, force = FALSE)

# GitHub
devtools::install_github("timothy-barry/ondisc")
devtools::install_github('satijalab/seurat-data')
devtools::install_git("git@github.com:Katsevich-Lab/lowmoi.git",
                      upgrade = "never")
