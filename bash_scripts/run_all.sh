#!/bin/sh

# 0. ensure system has required software dependencies:
#   - GCC version 11.1.0 or higher
#   - R version 4.1.2 or higher
#   - Python version 3.10.4 or higher
#   - Java version 8 or higher
#   - Nextflow version 21.10.6 or higher

###############
# PART 1: SETUP
###############
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

# 8. compute dataset sample sizes
Rscript ../R_scripts/compute_dataset_sample_sizes.R

#########################################################
# PART 2: RUN THE NEXTFLOW PIPELINES; PROCESS THE RESULTS
#########################################################
# 9. Undercover analysis
qsub ../pipeline_launch_scripts/undercover/grp_size_1_0523/grp_size_1.sh

# 10. Positive control analysis
qsub ../pipeline_launch_scripts/positive_control/pc_0523/pc_analysis.sh

# 11. MW resampling statistics analysis
qsub ../pipeline_launch_scripts/resampling_distributions/seurat_at_scale/seurat_at_scale.sh

# 12. Process the Nextflow pipeline results
Rscript ../R_scripts/process_results.R

####################################
# PART 3: RUN THE AUXILIARY ANALYSES 
####################################
# 13. run the analyses for figure s4
Rscript ../R_scripts/fig_s4_analyses.R

# 14. run the analyses for figure 5
Rscript ../R_scripts/save_datasets_as_r_objects.R
bash run_discovery_analyses.sh

# 15. run the analysis for supplementary figure s5
Rscript ../R_scripts/camp_simulation.R

# 16. compute the dataset statistical details
Rscript ../R_scripts/get_dataset_statistical_details.R

##########################
# PART 5: MAKE THE FIGURES
##########################
# fig 1
Rscript ../R_scripts/figure_creation/fig_1/fig_1.R
# fig 2
Rscript ../R_scripts/figure_creation/fig_2/fig_2.R # UPDATE!
# fig 3
Rscript ../R_scripts/figure_creation/fig_3/fig_3.R
# fig 4
Rscript ../R_scripts/figure_creation/fig_4/fig_4.R
# fig 5
Rscript ../R_scripts/figure_creation/fig_5/fig_5.R
# fig s1-s3
Rscript ../R_scripts/figure_creation/fig_s1_s3/fig_s1_s3.R
# fig s4
Rscript ../R_scripts/figure_creation/fig_s4/fig_s4.R
# fig s5
Rscript ../R_scripts/figure_creation/fig_s5/fig_s5.R