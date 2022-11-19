# limit NF driver to 5 GB of memory; source research config, pull current pipeline
export NXF_OPTS="-Xms500M -Xmx5G"
source ~/.research_config
nextflow pull https://github.com/Katsevich-Lab/undercover-grna-pipeline

# run pipeline
nextflow run Katsevich-Lab/undercover-grna-pipeline -r main \
 --data_method_pair_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/params_seurat_resampling.groovy" \
 --result_dir $LOCAL_SCEPTRE2_DATA_DIR"results/resampling_distributions" \
 --result_file_name "seurat_resampling_at_scale.rds" \
 --grna_modality "assignment" \
 --group_size "1" \
 --is_group_size_frac "false" \
 --partition_count "1" \
 --is_partition_count_frac "true" \
 --genes_to_subsample "2" \
 -profile aws \
 -with-trace \
 -resume

# is_partition_count_frac = false for trial
