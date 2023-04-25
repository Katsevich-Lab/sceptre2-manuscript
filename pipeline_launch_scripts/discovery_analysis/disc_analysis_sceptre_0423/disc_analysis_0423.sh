source ~/.research_config

nextflow pull Katsevich-Lab/pc-grna-pipeline
nextflow run pc-grna-pipeline -r main \
 --data_method_pair_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/discovery/params_sceptre_discovery_0423.groovy" \
 --grna_modality "assignment" \
 --result_dir $LOCAL_SCEPTRE2_DATA_DIR"/results/discovery_analyses" \
 --result_file_name "disc_analysis_sceptre_0423.rds" \
 --trial "false" \
 --pairs_file "tf_pairs_grouped.rds" \
 --grouped "true"