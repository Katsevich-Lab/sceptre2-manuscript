# CRAN
install.packages("readr", repos='http://cran.us.r-project.org')
install.packages("dplyr", repos='http://cran.us.r-project.org')
install.packages("readxl", repos='http://cran.us.r-project.org')
install.packages("Seurat", repos='http://cran.us.r-project.org')
install.packages("R.utils", repos='http://cran.us.r-project.org')
install.packages("purrr", repos='http://cran.us.r-project.org')
install.packages("RCurl", repos='http://cran.us.r-project.org')

# GitHub
devtools::install_github("Timothy-Barry/ondisc")
devtools::install_github('satijalab/seurat-data')
devtools::install_github("Katsevich-Lab/lowmoi")

# Bioconductor
BiocManager::install("GenomicRanges")

