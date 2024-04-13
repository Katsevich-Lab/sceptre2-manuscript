Reproducing the analyses reported in ‘Robust differential expression
testing for single-cell CRISPR screens at low multiplicty of infection’
================
Timothy Barry, Kaishu Mason, Kathryn Roeder, Eugene Katsevich

This repository contains code to reproduce the analyses reported in the
paper “Robust differential expression testing for single-cell CRISPR
screens at low multiplicty of infection” (Genome Biology, 2024). While
this repository serves as a useful starting point, given the complexity
of the codebase, we encourage readers interested in reproducing one or
more of the analyses to reach out to the authors (TB:
<tbarry@hsph.harvard.edu>, EK: <ekatsevi@wharton.upenn.edu>).

<!-- Data documentation is available [here](https://github.com/Katsevich-Lab/sceptre2-manuscript/blob/main/docs/data_doc.pdf). -->

# Dependencies

Ensure that your system has the required software dependencies.

- GCC version 11.1.0 or higher
- R version 4.1.2 or higher
- Python version 3.10.4 or higher
- Java version 8 or higher
- Nextflow version 21.10.6 or higher

# Get started

First, clone the `sceptre2-manuscript` repository onto your machine.

    git clone git@github.com:Katsevich-Lab/sceptre2-manuscript.git

Next, manually download the processed datasets from Dropbox. (You
alternately can automatically download and process the datasets via the
command `bash import_all_data.sh`, but this is slow and potentially
could hit a few snags depending on your system.) Download the following
three directories from the [Dropbox data
repository](www.dropbox.com/sh/jekmk1v4mr4kj3b/AAAhznGqk-TIZKhW40xiU6ORa?dl=0):
`frangieh-2021`, `papalexi-2021`, `schraivogel-2020`. The code
originally used to download and process these datasets (alongside
helpful documentation) is available in the following repositories:
[Frangieh](https://github.com/Katsevich-Lab/import-frangieh-2021),
[Papalexi](https://github.com/Katsevich-Lab/import-papalexi-2021),
[Schraivogel](https://github.com/Katsevich-Lab/import-schraivogel-2020).

We used a config file to increase the portability of our code across
machines. Create a config file called `.research_config` in your home
directory.

    cd
    touch .research_config

Define the following variables within this file:

- `LOCAL_PAPALEXI_2021_DATA_DIR`: the location of the `papalexi-2021`
  data directory

- `LOCAL_SCHRAIVOGEL_2020_DATA_DIR`: the location of the
  `schraivogel-2020` data directory

- `LOCAL_FRANGIEH_2021_DATA_DIR`: the location of the `frangieh-2021`
  data directory

- `LOCAL_CODE_DIR`: the location of the directory in which the cloned
  `sceptre2-manuscript` is located

- `LOCAL_SCEPTRE2_DATA_DIR`: the location of the directory in which to
  store files (e.g., auxiliary data, results) associated with the paper.

- `KATS_SOFT_DIR`: the location of the directory in which to create the
  Python virtual environment

The contents of the `.research_config` file should look like something
along the following lines.

    LOCAL_PAPALEXI_2021_DATA_DIR=/Users/timbarry/research_offsite/external/papalexi-2021/
    LOCAL_SCHRAIVOGEL_2020_DATA_DIR=/Users/timbarry/research_offsite/external/schraivogel-2020/
    LOCAL_FRANGIEH_2021_DATA_DIR=/Users/timbarry/research_offsite/external/frangieh-2021/
    LOCAL_CODE_DIR=/Users/timbarry/research_code/
    LOCAL_SCEPTRE2_DATA_DIR=/Users/timbarry/research_offsite/projects/sceptre2/
    KATS_SOFT_DIR=/Users/timbarry/

Next, create an `.Rprofile` file in your home directory (if you have not
yet done so).

    cd
    touch .Rprofile

Add the following commands to your `.Rprofile`.

    suppressMessages(library(conflicted))
    .get_config_path <- function(dir_name) {
      cmd <- paste0("source ~/.research_config; echo $", dir_name)
      system(command = cmd, intern = TRUE)
    }
    Sys.setenv(RETICULATE_PYTHON = paste0(.get_config_path("KATS_SOFT_DIR"), "py/lowmoi-venv/bin/python"))
    utils::setRepositories(ind = 1:4)

# Set up the analyses

Navigate to the `bash_scripts` subdirectory of the `sceptre2-manuscript`
directory. Execute the following command to create a Python virtual
environment called `py/lowmoi-venv` and download the necessary Python
packages.

    cd bash_scripts
    bash install_python_packages.sh

If you are having trouble initializing the Python virtual environment,
consult the [following
documentation](https://github.com/Katsevich-Lab/sceptre2-manuscript/blob/main/docs/setting_up_python.pdf).
Next, install the required R packages.

    Rscript ../R_scripts/install_R_packages.R

Finally, complete the rest of the setup, which consists of initializing
the directory structure of the SCEPTRE2 data directory, initializing the
symbolic links, creating the synthetic data, performing quality control
on the data, and computing sample sizes. The following code takes about
20 minutes to execute.

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

We recommend downloading the results from the [Dropbox results
repository](https://www.dropbox.com/sh/76sviudz0cg0mgy/AACoigzD0CWGa9S5HVgG2gHua?dl=0),
so that you can reproduce the figures without having to rerun the
analyses. This can be executed via the following commands.

    source ~/.research_config
    cd $LOCAL_SCEPTRE2_DATA_DIR"/results"
     wget --max-redirect=20 -O download.zip https://www.dropbox.com/sh/yuubaoro3k75c61/AACPi9LDjY0B1pxwMxUMSoTca?dl=1
     unzip -o download.zip

All results will be placed into the correct locations.

# Run the Nextflow pipelines

The next step is to run a sequence of Nextflow pipelines to carry out
the main analyses. These pipelines should be run one-by-one.

    # 9. Undercover analysis
    qsub ../pipeline_launch_scripts/undercover/grp_size_1_0523/grp_size_1.sh

    # 10. Positive control analysis
    qsub ../pipeline_launch_scripts/positive_control/pc_0523/pc_analysis.sh

    # 11. MW resampling statistics analysis
    qsub ../pipeline_launch_scripts/resampling_distributions/seurat_at_scale/seurat_at_scale.sh

    # 12. Unfiltered SCEPTRE positive control analysis
    qsub ../pipeline_launch_scripts/positive_control/sceptre_unfiltered_pc_0523/pc_analysis.sh

We carried out our analyses on an SGE cluster; thus, we submitted the
Nextflow jobs to the scheduler via `qsub`. If you instead are running
the analyses on a SLURM cluster, for example, then you should submit the
jobs via `sbatch`. If you are running the analyses locally, then you
should submit the jobs via `bash`, etc.

Once the jobs have completed, process the results.

    Rscript ../R_scripts/process_results.R

# Run the auxiliary analyses

The next step is to run several auxiliary analyses. These analyses are
fairly small and thus do not necessitate Nextflow pipelines.

    # 14. run the analyses for figure s4
    Rscript ../R_scripts/fig_s4_analyses.R

    # 15. run the analyses for figure 5
    Rscript ../R_scripts/save_datasets_as_r_objects.R
    bash run_discovery_analyses.sh

    # 16. run the analysis for supplementary figure s5
    Rscript ../R_scripts/camp_simulation.R

    # 17. compute the dataset statistical details
    Rscript ../R_scripts/get_dataset_statistical_details.R

# Create the figures

The final step is to create the figures. To do so issue the following
commands.

    # fig 1
    Rscript ../R_scripts/figure_creation/fig_1/fig_1.R
    # fig 2
    Rscript ../R_scripts/figure_creation/fig_2/fig_2.R
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
    # fig s6
    Rscript ../R_scripts/figure_creation/fig_s6/fig_s6.R
    # fig s7
    Rscript ../R_scripts/figure_creation/fig_s7/fig_s7.R

The script `bash/run_all.sh` runs the entire analysis from start to
finish, from downloading and processing the raw data to creating the
figures.

Contact the authors if you seek help with reproducing an aspect of the
analysis.

- Timothy (Tim) Barry: tbarry2@andrew.cmu.edu

- Eugene (Gene) Katsevich: ekatsevi@wharton.upenn.edu
