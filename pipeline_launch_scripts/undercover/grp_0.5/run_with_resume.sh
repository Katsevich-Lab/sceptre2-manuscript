curr_time=$(date '+%m%d%H%M')

nextflow run undercover-grna-pipeline -r main \
 --data_method_pair_file /home/stat/timbar/research_code/sceptre2-manuscript/param_files/data_method_pair_file.groovy \
 --result_dir /home/stat/timbar/data/projects/sceptre2/results/undercover_grna_analysis \
 --result_file_name undercover_result_grp_size_half.rds \
 --grna_modality assignment \
 --group_size 0.5 \
 --is_group_size_frac true \
 --partition_count 1 \
 --is_partition_count_frac true \
 --machine_name hpcc \
 --time $curr_time \
 -profile standard
