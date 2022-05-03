source ~/.research_config

# Initialize several offsite directories within sceptre2
# |- data
#   |- schraivogel
#      |- tap
#      |- perturb
#   |- papalexi
#      |- gene
#      |- protein
#   |- simulated
#      |- gene
#      |- gRNA
# |- results

sceptre2_dir=$LOCAL_SCEPTRE2_DATA_DIR
mkdir -p $sceptre2_dir"data/simulated" \
 $sceptre2_dir"data/schraivogel/tap" \
 $sceptre2_dir"data/schraivogel/perturb" \
 $sceptre2_dir"data/papalexi/gene" \
 $sceptre2_dir"data/papalexi/protein" \
 $sceptre2_dir"data/simulated/gene" \
 $sceptre2_dir"data/simulated/gRNA" \
 $sceptre2_dir"results"
