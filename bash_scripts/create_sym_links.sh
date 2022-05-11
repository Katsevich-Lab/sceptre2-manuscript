# creates sym links to backing .odm files

source ~/.research_config
sceptre2_data_dir=$LOCAL_SCEPTRE2_DATA_DIR"data/"
papalexi_proc_dir=$LOCAL_PAPALEXI_2021_DATA_DIR"processed/"
liscovitch_proc_dir=$LOCAL_LISCOVITCH_2021_DATA_DIR"processed/"
schraivogel_proc_dir=$LOCAL_SCHRAIVOGEL_2020_DATA_DIR"processed/"

# Papalexi data
# gene
ln -s $papalexi_proc_dir"gene/expression_matrix.odm" $sceptre2_data_dir"papalexi/eccite_screen/gene/matrix.odm"
ln -s $papalexi_proc_dir"gene/metadata.rds" $sceptre2_data_dir"papalexi/eccite_screen/gene/metadata_orig.rds"
# grna
ln -s $papalexi_proc_dir"gRNA/count_matrix.odm" $sceptre2_data_dir"papalexi/eccite_screen/grna/matrix.odm"
ln -s $papalexi_proc_dir"gRNA/metadata.rds" $sceptre2_data_dir"papalexi/eccite_screen/grna/metadata_orig.rds"
# protein
ln -s $papalexi_proc_dir"protein/count_matrix.odm" $sceptre2_data_dir"papalexi/eccite_screen/protein/matrix.odm"
ln -s $papalexi_proc_dir"protein/metadata.rds" $sceptre2_data_dir"papalexi/eccite_screen/protein/metadata_orig.rds"

# Liscovitch data
# experiment_small chromatin
ln -s $liscovitch_proc_dir"experiment_small/chromatin/chip_counts.odm" $sceptre2_data_dir"liscovitch/experiment_small/chromatin/matrix.odm"
ln -s $liscovitch_proc_dir"experiment_small/chromatin/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_small/chromatin/metadata_orig.rds"
# experiment_small gRNA
ln -s $liscovitch_proc_dir"experiment_small/gRNA/gRNA_counts.odm" $sceptre2_data_dir"liscovitch/experiment_small/gRNA/matrix.odm"
ln -s $liscovitch_proc_dir"experiment_small/gRNA/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_small/gRNA/metadata_orig.rds"
# experiment_big chromatin
ln -s $liscovitch_proc_dir"experiment_big/chromatin/chip_counts.odm" $sceptre2_data_dir"liscovitch/experiment_big/chromatin/matrix.odm"
ln -s $liscovitch_proc_dir"experiment_big/chromatin/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_big/chromatin/metadata_orig.rds"
# experiment_big gRNA
ln -s $liscovitch_proc_dir"experiment_big/gRNA/gRNA_counts.odm" $sceptre2_data_dir"liscovitch/experiment_big/gRNA/matrix.odm"
ln -s $liscovitch_proc_dir"experiment_big/gRNA/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_big/gRNA/metadata_orig.rds"

# Schraivogel data
# enhancer_screen_chr11 gene
ln -s $schraivogel_proc_dir"enhancer_screen_chr11/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/gene/matrix.odm"
ln -s $schraivogel_proc_dir"enhancer_screen_chr11/gene/metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/gene/metadata_orig.rds"
# enhancer_screen_chr11 gRNA
ln -s $schraivogel_proc_dir"enhancer_screen_chr11/gRNA/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/gRNA/matrix.odm"
ln -s $schraivogel_proc_dir"enhancer_screen_chr11/gRNA/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/gRNA/metadata_orig.rds"
# enhancer_screen_chr8 gene
ln -s $schraivogel_proc_dir"enhancer_screen_chr8/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/gene/matrix.odm"
ln -s $schraivogel_proc_dir"enhancer_screen_chr8/gene/metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/gene/metadata_orig.rds"
# enhancer_screen_chr8 gRNA
ln -s $schraivogel_proc_dir"enhancer_screen_chr8/gRNA/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/gRNA/matrix.odm"
ln -s $schraivogel_proc_dir"enhancer_screen_chr8/gRNA/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/gRNA/metadata_orig.rds"
# ground_truth_perturbseq gene
ln -s $schraivogel_proc_dir"ground_truth_perturbseq/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/gene/matrix.odm"
ln -s $schraivogel_proc_dir"ground_truth_perturbseq/gene/metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/gene/metadata_orig.rds"
# ground_truth_perturbseq gRNA
ln -s $schraivogel_proc_dir"ground_truth_perturbseq/gRNA/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/gRNA/matrix.odm"
ln -s $schraivogel_proc_dir"ground_truth_perturbseq/gRNA/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/gRNA/metadata_orig.rds"
# ground_truth_tapseq gene
ln -s $schraivogel_proc_dir"ground_truth_tapseq/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/gene/matrix.odm"
ln -s $schraivogel_proc_dir"ground_truth_tapseq/gene/metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/gene/metadata_orig.rds"
# ground_truth_tapseq gRNA
ln -s $schraivogel_proc_dir"ground_truth_tapseq/gRNA/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/gRNA/matrix.odm"
ln -s $schraivogel_proc_dir"ground_truth_tapseq/gRNA/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/gRNA/metadata_orig.rds"

# Franghei
