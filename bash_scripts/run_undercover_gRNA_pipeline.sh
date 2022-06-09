#$ -j y

export NXF_OPTS="-Xms500M -Xmx3G" # limit NF driver to 3 GB of memory
source ~/.research_config

curr_time=$(date '+%m%d%H%M')
#nextflow pull https://github.com/Katsevich-Lab/undercover-gRNA-pipeline
#nextflow run undercover-gRNA-pipeline -r main \
nextflow run $LOCAL_CODE_DIR"/undercover-gRNA-pipeline/main.nf" \
 --data_method_pair_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/test_file.groovy" \
 --result_dir $LOCAL_SCEPTRE2_DATA_DIR"results/undercover_grna_analysis" \
 --result_file_name "trial.rds" \
 --one_neg_control "TRUE" \
 --max_retries 4 \
 --grna_modality "assignment" \
 --machine_name $MACHINE_NAME \
 --time $curr_time
