source ~/.research_config

nextflow pull Katsevich-Lab/pc-grna-pipeline
nextflow run pc-grna-pipeline -r main \
 --data_method_pair_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/pc/params_pc_seuratnb_0423.groovy" \
 --grna_modality "assignment" \
 --result_dir $LOCAL_SCEPTRE2_DATA_DIR"/results/positive_control_analysis" \
 --result_file_name "pc_results_seuratnb_0423.rds" \
 --trial "false" \
 --pairs_file "pos_control_pairs_grouped.rds" \
 --grouped "true"