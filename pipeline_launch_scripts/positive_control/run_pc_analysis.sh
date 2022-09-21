source ~/.research_config

nextflow pull Katsevich-Lab/pc-grna-pipeline
nextflow run pc-grna-pipeline -r main \
 --data_method_pair_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/data_method_pairs_pc.groovy" \
 --grna_modality "assignment" \
 --result_dir $LOCAL_SCEPTRE2_DATA_DIR"/results/positive_control_analysis" \
 --result_file_name "pc_results.rds" \
 --trial "true" \
 -resume
