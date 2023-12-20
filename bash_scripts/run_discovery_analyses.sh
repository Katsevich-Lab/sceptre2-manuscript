source ~/.research_config
result_dir=$LOCAL_SCEPTRE2_DATA_DIR"results/discovery_analyses/"

for dataset in papalexi frangieh; do
  for full_test_statistic in FALSE TRUE; do
    out_file=$result_dir$dataset"_full_stat_"$full_test_statistic"_benchmark.out"
    /usr/bin/time -l -h -p ../R_scripts/run_discovery_analysis.R $dataset $full_test_statistic 2> $out_file
  done;
done
