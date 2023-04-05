source ~/.research_config
result_dir=$LOCAL_SCEPTRE2_DATA_DIR"results/discovery_analyses/"

for dataset in papalexi frangieh; do
  for analysis_type in calibration discovery; do
    out_file=$result_dir$dataset"_"$analysis_type"_benchmark.out"
    /usr/bin/time -l -h -p ../R_scripts/run_discovery_analysis.R $dataset $analysis_type 2> $out_file
  done;
done