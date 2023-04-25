# limit NF driver to 5 GB of memory; source research config, pull current pipeline
export NXF_OPTS="-Xms500M -Xmx5G"
source ~/.research_config
nextflow pull https://github.com/Katsevich-Lab/undercover-grna-pipeline

# run pipeline
nextflow run Katsevich-Lab/undercover-grna-pipeline -r main \
 --data_method_pair_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/undercover/params_undercover_seuratnb_0423.groovy" \
 --result_dir $LOCAL_SCEPTRE2_DATA_DIR"results/undercover_grna_analysis" \
 --result_file_name "undercover_result_grp_1_seuratnb_0423.rds" \
 --grna_modality "assignment" \
 --group_size "1" \
 --is_group_size_frac "false" \
 --partition_count "1" \
 --is_partition_count_frac "true" \
 -profile standard \
 -with-trace

# is_partition_count_frac = false for trial
