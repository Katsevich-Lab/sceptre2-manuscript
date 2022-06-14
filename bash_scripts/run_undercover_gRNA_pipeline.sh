#$ -j y

export NXF_OPTS="-Xms500M -Xmx4G" # limit NF driver to 4 GB of memory
source ~/.research_config

curr_time=$(date '+%m%d%H%M')
nextflow pull https://github.com/Katsevich-Lab/undercover-gRNA-pipeline
nextflow run undercover-gRNA-pipeline -r main \
 --data_method_pair_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/data_method_pair_file.groovy" \
 --result_dir $LOCAL_SCEPTRE2_DATA_DIR"results/undercover_grna_analysis" \
 --result_file_name "result_mimosca.rds" \
 --one_neg_control "FALSE" \
 --max_retries 2 \
 --grna_modality "assignment" \
 --machine_name $MACHINE_NAME \
 --time $curr_time \
 -with-trace
