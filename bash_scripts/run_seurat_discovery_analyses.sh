source ~/.research_config

datasets=(schraivogel/enhancer_screen_chr11/gene schraivogel/enhancer_screen_chr8/gene papalexi/eccite_screen/gene frangieh/control/gene)
for dataset in ${datasets[@]}; do
  qsub -j y -l m_mem_free=10G run_seurat_discovery_analysis.sh $dataset
done
