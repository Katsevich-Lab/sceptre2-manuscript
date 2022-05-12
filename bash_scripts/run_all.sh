#!/bin/sh

# 1. set up the offsite directory structure
bash setup.sh
# 2. create the symbolic links within the offsite directory structure
bash create_sym_links.sh
# 3. cell-specific QC on Frangieh and Liscovitch data
Rscript ../R_scripts/cell_qc.R
# 4. feature-specific QC on all data
Rscript ../R_scripts/feature_qc.R
# 5. create the synthetic data
Rscript ../R_scripts/create_synthetic_data.R
