source ~/.research_config

# Each leaf directory of the data directory has the following form:
# data/paper_name/dataset_name/modality/metadata_qc.rds
# paper_name: one of liscovitch, schraivogel, papalexi, frangieh
# dataset_name: a descriptive name of the dataset
# modality: the modality within the dataset; one of gene, grna, protein, and chromatin

# Initialize several offsite directories within sceptre2
# |- data
#   |- schraivogel
#      |- ground_truth_tapseq
#          |- gene
#          |- grna
#      |- ground_truth_perturbseq
#          |- gene
#          |- grna
#      |- enhancer_screen_chr8
#          |- gene
#          |- grna
#      |- enhancer_screen_chr11
#          |- gene
#          |- grna
#   |- papalexi
#      |- eccite_screen
#          |- gene
#          |- grna
#          |- protein
#   |- liscovitch
#       |- experiment_small
#          |- chromatin
#          |- grna
#       |- experiment_big
#          |- chromatin
#          |- grna
#   |- frangieh
#       |- control
#          |- gene
#          |- protein
#          |- grna
#       |- ifn-gamma
#          |- gene
#          |- protein
#          |- grna
#       |- co-culture
#          |- gene
#          |- protein
#          |- grna
#   |- simulated
#      |- gene
#      |- grna
# |- results

sceptre2_dir=$LOCAL_SCEPTRE2_DATA_DIR
sceptre2_data_dir=$sceptre2_dir"data/"

rm -rf $sceptre2_data_dir
mkdir -p $sceptre2_data_dir"schraivogel/ground_truth_tapseq/gene" \
 $sceptre2_data_dir"schraivogel/ground_truth_tapseq/grna" \
 $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/gene" \
 $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/grna" \
 $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/gene" \
 $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/grna" \
 $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/gene" \
 $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/grna" \
 $sceptre2_data_dir"papalexi/eccite_screen/gene" \
 $sceptre2_data_dir"papalexi/eccite_screen/grna" \
 $sceptre2_data_dir"papalexi/eccite_screen/protein" \
 $sceptre2_data_dir"liscovitch/experiment_small/chromatin" \
 $sceptre2_data_dir"liscovitch/experiment_small/grna" \
 $sceptre2_data_dir"liscovitch/experiment_big/grna" \
 $sceptre2_data_dir"liscovitch/experiment_big/chromatin" \
 $sceptre2_data_dir"frangieh/co_culture/gene" \
 $sceptre2_data_dir"frangieh/co_culture/protein" \
 $sceptre2_data_dir"frangieh/co_culture/grna" \
 $sceptre2_data_dir"frangieh/control/gene" \
 $sceptre2_data_dir"frangieh/control/protein" \
 $sceptre2_data_dir"frangieh/control/grna" \
 $sceptre2_data_dir"frangieh/ifn_gamma/gene" \
 $sceptre2_data_dir"frangieh/ifn_gamma/protein" \
 $sceptre2_data_dir"frangieh/ifn_gamma/grna" \
 $sceptre2_data_dir"simulated/experiment_1/gene" \
 $sceptre2_data_dir"simulated/experiment_1/grna" \
 $sceptre2_dir"results"
