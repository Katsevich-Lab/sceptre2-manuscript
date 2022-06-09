#!/bin/sh

# 0. ensure system has required software dependencies:
#   - GCC version 11.1.0 or higher
#   - R version 4.1.2 or higher
#   - Python version 3.10.4 or higher
#   - Java version 8 or higher
#   - Nextflow version 21.10.6 or higher

# 1. install required Python packages
bash install_python_packages.sh

# 2. install required R packages
Rscript ../R_scripts/install_R_packages.R

# 3. import all datasets
bash import_all_data.sh

# 4. set up the offsite directory structure
bash setup.sh

# 5. create the symbolic links within the offsite directory structure
bash create_sym_links.sh

# 6. create the synthetic data
Rscript ../R_scripts/create_synthetic_data.R

# 7. run the uniform processing step (which includes cell and feature qc)
Rscript ../R_scripts/uniform_processing.R
