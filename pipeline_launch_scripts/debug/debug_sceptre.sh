# limit NF driver to 5 GB of memory; source research config, pull current pipeline
export NXF_OPTS="-Xms500M -Xmx5G"
source ~/.research_config
nextflow pull https://github.com/Katsevich-Lab/undercover-grna-pipeline
curr_time=$(date '+%m%d%H%M')

# run pipeline
nextflow run Katsevich-Lab/undercover-grna-pipeline -r main \
 --data_method_pair_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/sceptre_debug.groovy" \
 --result_dir $LOCAL_SCEPTRE2_DATA_DIR"results/undercover_grna_analysis" \
 --result_file_name "sceptre_debug_kde.rds" \
 --grna_modality "assignment" \
 --group_size "1" \
 --is_group_size_frac "false" \
 --partition_count 1 \
 --is_partition_count_frac "true" \
 --machine_name $MACHINE_NAME \
 --time $curr_time \
 -profile standard
