#!/bin/bash

# Limit NF driver to 5 GB memory
export NXF_OPTS="-Xms500M -Xmx5G"

#############
# INPUT FILES
#############
source ~/.research_config
gasp_data_dir=$LOCAL_SCEPTRE2_DATA_DIR"data/gasperini/at_scale/"
# i) multimodal metadata file
multimodal_metadata_fp=$gasp_data_dir"multimodal_metadata.rds"
# ii) gene ODM
gene_odm_fp=$gasp_data_dir"gene/matrix.odm"
# iii) grna ODM
grna_odm_fp=$gasp_data_dir"grna_expression/matrix.odm"
# iv) gene-grna group pairs
pair_fp=$LOCAL_SCEPTRE2_DATA_DIR/data/gasperini/at_scale/benchmarking_pairs.rds

if [ ! -f "$pair_fp" ]; then
  echo "Generating benchmark pairs..."
  Rscript $LOCAL_CODE_DIR/sceptre2-manuscript/writeups/resampling_distributions/create_benchmark_pairs.R
fi

###############
# OPTIONAL ARGS
###############
# formula, threshold, B, side, n_pairs_to_sample, gene_pod_size, grna_group_pod_size, pair_pod_size are optional args
formula="~p_mito+batch+log(gene_n_nonzero)+log(gene_n_umis)+log(grna_expression_n_nonzero)+log(grna_expression_n_umis)"
gene_pod_size=400
grna_group_pod_size=400
pair_pod_size=400
# gene_pod_size=4
# grna_group_pod_size=4
# pair_pod_size=4

########################
# invoke the NF pipeline
########################

nextflow pull ekatsevi/sceptre-pipeline
Rscript -e 'devtools::install_github("ekatsevi/sceptre")'

for method in "gcm" "crt" "gcm_crt"; do
  result_fp=$LOCAL_SCEPTRE2_DATA_DIR/results/highmoi_pipeline/gasperini_$method"_result.rds"
  if [ ! -f $result_fp ]; then
    echo "Running "$method" analysis..."
    nextflow run ekatsevi/sceptre-pipeline -r main \
     --multimodal_metadata_fp $multimodal_metadata_fp \
     --gene_odm_fp $gene_odm_fp \
     --grna_odm_fp $grna_odm_fp \
     --pair_fp $pair_fp \
     --result_fp $result_fp \
     --formula $formula \
     --grna_modality_name "grna_expression" \
     --inference_method $method \
     --gene_pod_size $gene_pod_size \
     --grna_group_pod_size $grna_group_pod_size \
     --pair_pod_size $pair_pod_size \
     -profile local \
     -resume
  else
    echo $method" analysis already done!"
  fi
done 