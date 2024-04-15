source ~/.research_config
result_dir=$LOCAL_SCEPTRE2_DATA_DIR"results/discovery_analyses/fig_s12/"

for dataset in papalexi frangieh; do
  for full_test_statistic in FALSE TRUE; do
    ../R_scripts/run_discovery_analysis_fig_s12.R $dataset $full_test_statistic
  done;
done
