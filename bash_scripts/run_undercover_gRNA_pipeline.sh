export NXF_OPTS="-Xms500M -Xmx3G" # limit NF to 3 GB of memory
source ~/.research_config

curr_time=$(date '+%m%d%H%M')
nextflow pull https://github.com/Katsevich-Lab/undercover-gRNA-pipeline
nextflow run undercover-gRNA-pipeline -r main \
 --data_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/data_file.R" \
 --data_method_pair_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/data_method_pair_file.groovy" \
 --machine_name $MACHINE_NAME \
 --result_dir $LOCAL_SCEPTRE2_DATA_DIR"results" \
 --time $curr_time
