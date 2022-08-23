#$ -j y

export NXF_OPTS="-Xms500M -Xmx5G" # limit NF driver to 5 GB of memory
source ~/.research_config
nextflow pull https://github.com/Katsevich-Lab/undercover-grna-pipeline

# group_sizes=("1" "2" "3" "0.5")
# is_group_size_fracs=("false" "false" "false" "true")
# result_file_names=("undercover_result_grp_size_1.rds" "undercover_result_grp_size_2.rds" "undercover_result_grp_size_3.rds" "undercover_result_grp_size_half.rds")

group_sizes=("1")
is_group_size_fracs=("false")
result_file_names=("undercover_result_grp_size_1_sceptre.rds")

for i in ${!group_sizes[@]}; do
  cur_dir="grp_${group_sizes[$i]}"
  mkdir -p $cur_dir
  cd $cur_dir

  curr_time=$(date '+%m%d%H%M')
  nextflow run undercover-grna-pipeline -r main \
   --data_method_pair_file $LOCAL_CODE_DIR"/sceptre2-manuscript/param_files/data_method_pair_file.groovy" \
   --result_dir $LOCAL_SCEPTRE2_DATA_DIR"results/undercover_grna_analysis" \
   --result_file_name "${result_file_names[$i]}" \
   --grna_modality "assignment" \
   --group_size "${group_sizes[$i]}" \
   --is_group_size_frac "${is_group_size_fracs[$i]}" \
   --partition_count 1 \
   --is_partition_count_frac "true" \
   --machine_name $MACHINE_NAME \
   --time $curr_time \
   -profile standard
  wait
  cd ..
done
