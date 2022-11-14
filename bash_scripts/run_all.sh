#!/bin/sh

# 0. ensure system has required software dependencies:
#   - GCC version 11.1.0 or higher
#   - R version 4.1.2 or higher
#   - Python version 3.10.4 or higher
#   - Java version 8 or higher
#   - Nextflow version 21.10.6 or higher

#######################################
# PART 1: SETUP
#######################################

# 1. install required Python packages
bash install_python_packages.sh

# 2. install required R packages
Rscript ../R_scripts/install_R_packages.R

# 3. import all datasets
bash import_all_data.sh
wait

# 4. set up the offsite directory structure
bash setup.sh

# 5. create the symbolic links within the offsite directory structure
bash create_sym_links.sh

# 6. create the synthetic data
Rscript ../R_scripts/create_synthetic_data.R

# 7. run the low MOI uniform processing step (which includes cell and feature qc)
Rscript ../R_scripts/uniform_processing_lowmoi.R

# 8. run the high MOI uniform processing step
Rscript ../R_scripts/uniform_processing_highmoi.R

#######################################
# PART 2: RESULTS GENERATION
#######################################
# 9. Undercover results for group size = 1, all methods and datasets
qsub ../pipeline_launch_scripts/undercover/grp_size_1/grp_size_1.sh

# 10. Undercover results for group size = 2, all methods and datasets

# 11. Undercover results for group size = half, all methods and datasets


#######################################
# PART 3: MAKING FIGURES
#######################################
