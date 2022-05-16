# creates sym links to backing .odm files

source ~/.research_config
sceptre2_data_dir=$LOCAL_SCEPTRE2_DATA_DIR"data/"
papalexi_proc_dir=$LOCAL_PAPALEXI_2021_DATA_DIR"processed/"
liscovitch_proc_dir=$LOCAL_LISCOVITCH_2021_DATA_DIR"processed/"
schraivogel_proc_dir=$LOCAL_SCHRAIVOGEL_2020_DATA_DIR"processed/"
frangieh_proc_dir=$LOCAL_FRANGIEH_2021_DATA_DIR"processed/perturb-cite-seq/"

# Papalexi data
# gene
ln -s $papalexi_proc_dir"gene/expression_matrix.odm" $sceptre2_data_dir"papalexi/eccite_screen/gene/matrix.odm"
ln -s $papalexi_proc_dir"gene/metadata.rds" $sceptre2_data_dir"papalexi/eccite_screen/gene/metadata_orig.rds"
# grna
ln -s $papalexi_proc_dir"grna/count_matrix.odm" $sceptre2_data_dir"papalexi/eccite_screen/grna/matrix.odm"
ln -s $papalexi_proc_dir"grna/metadata.rds" $sceptre2_data_dir"papalexi/eccite_screen/grna/metadata_orig.rds"
# protein
ln -s $papalexi_proc_dir"protein/count_matrix.odm" $sceptre2_data_dir"papalexi/eccite_screen/protein/matrix.odm"
ln -s $papalexi_proc_dir"protein/metadata.rds" $sceptre2_data_dir"papalexi/eccite_screen/protein/metadata_orig.rds"

# Liscovitch data
# experiment_small chromatin
ln -s $liscovitch_proc_dir"experiment_small/chromatin/chip_counts.odm" $sceptre2_data_dir"liscovitch/experiment_small/chromatin/matrix.odm"
ln -s $liscovitch_proc_dir"experiment_small/chromatin/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_small/chromatin/metadata_orig.rds"
# experiment_small grna
ln -s $liscovitch_proc_dir"experiment_small/grna/grna_counts.odm" $sceptre2_data_dir"liscovitch/experiment_small/grna/matrix.odm"
ln -s $liscovitch_proc_dir"experiment_small/grna/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_small/grna/metadata_orig.rds"
# experiment_big chromatin
ln -s $liscovitch_proc_dir"experiment_big/chromatin/chip_counts.odm" $sceptre2_data_dir"liscovitch/experiment_big/chromatin/matrix.odm"
ln -s $liscovitch_proc_dir"experiment_big/chromatin/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_big/chromatin/metadata_orig.rds"
# experiment_big grna
ln -s $liscovitch_proc_dir"experiment_big/grna/grna_counts.odm" $sceptre2_data_dir"liscovitch/experiment_big/grna/matrix.odm"
ln -s $liscovitch_proc_dir"experiment_big/grna/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_big/grna/metadata_orig.rds"

# Schraivogel data
# enhancer_screen_chr11 gene
ln -s $schraivogel_proc_dir"enhancer_screen_chr11/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/gene/matrix.odm"
ln -s $schraivogel_proc_dir"enhancer_screen_chr11/gene/metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/gene/metadata_orig.rds"
# enhancer_screen_chr11 grna
ln -s $schraivogel_proc_dir"enhancer_screen_chr11/grna/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/grna/matrix.odm"
ln -s $schraivogel_proc_dir"enhancer_screen_chr11/grna/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/grna/metadata_orig.rds"
# enhancer_screen_chr8 gene
ln -s $schraivogel_proc_dir"enhancer_screen_chr8/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/gene/matrix.odm"
ln -s $schraivogel_proc_dir"enhancer_screen_chr8/gene/metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/gene/metadata_orig.rds"
# enhancer_screen_chr8 grna
ln -s $schraivogel_proc_dir"enhancer_screen_chr8/grna/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/grna/matrix.odm"
ln -s $schraivogel_proc_dir"enhancer_screen_chr8/grna/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/grna/metadata_orig.rds"
# ground_truth_perturbseq gene
ln -s $schraivogel_proc_dir"ground_truth_perturbseq/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/gene/matrix.odm"
ln -s $schraivogel_proc_dir"ground_truth_perturbseq/gene/metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/gene/metadata_orig.rds"
# ground_truth_perturbseq grna
ln -s $schraivogel_proc_dir"ground_truth_perturbseq/grna/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/grna/matrix.odm"
ln -s $schraivogel_proc_dir"ground_truth_perturbseq/grna/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/grna/metadata_orig.rds"
# ground_truth_tapseq gene
ln -s $schraivogel_proc_dir"ground_truth_tapseq/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/gene/matrix.odm"
ln -s $schraivogel_proc_dir"ground_truth_tapseq/gene/metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/gene/metadata_orig.rds"
# ground_truth_tapseq grna
ln -s $schraivogel_proc_dir"ground_truth_tapseq/grna/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/grna/matrix.odm"
ln -s $schraivogel_proc_dir"ground_truth_tapseq/grna/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/grna/metadata_orig.rds"

# Franghei
# copy the gene odms
ln -s $frangieh_proc_dir"gene/gene_expression_matrix.odm" $sceptre2_data_dir"frangieh/co_culture/gene/matrix.odm"
ln -s $frangieh_proc_dir"gene/gene_expression_matrix.odm" $sceptre2_data_dir"frangieh/control/gene/matrix.odm"
ln -s $frangieh_proc_dir"gene/gene_expression_matrix.odm" $sceptre2_data_dir"frangieh/ifn_gamma/gene/matrix.odm"
# copy the grna odms
ln -s $frangieh_proc_dir"grna/grna_assignments_ungrouped.odm" $sceptre2_data_dir"frangieh/co_culture/grna/matrix.odm"
ln -s $frangieh_proc_dir"grna/grna_assignments_ungrouped.odm" $sceptre2_data_dir"frangieh/control/grna/matrix.odm"
ln -s $frangieh_proc_dir"grna/grna_assignments_ungrouped.odm" $sceptre2_data_dir"frangieh/ifn_gamma/grna/matrix.odm"
# copy the protein odms
ln -s $frangieh_proc_dir"protein/protein_expression_matrix.odm" $sceptre2_data_dir"frangieh/co_culture/protein/matrix.odm"
ln -s $frangieh_proc_dir"protein/protein_expression_matrix.odm" $sceptre2_data_dir"frangieh/control/protein/matrix.odm"
ln -s $frangieh_proc_dir"protein/protein_expression_matrix.odm" $sceptre2_data_dir"frangieh/ifn_gamma/protein/matrix.odm"
# copy the gene metadata
ln -s $frangieh_proc_dir"gene/gene_expression_metadata_co-culture.rds" $sceptre2_data_dir"frangieh/co_culture/gene/metadata_orig.rds"
ln -s $frangieh_proc_dir"gene/gene_expression_metadata_control.rds" $sceptre2_data_dir"frangieh/control/gene/metadata_orig.rds"
ln -s $frangieh_proc_dir"gene/gene_expression_metadata_ifn-gamma.rds" $sceptre2_data_dir"frangieh/ifn_gamma/gene/metadata_orig.rds"
# copy the grna metadata
ln -s $frangieh_proc_dir"grna/grna_assignments_ungrouped_metadata_co-culture.rds" $sceptre2_data_dir"frangieh/co_culture/grna/metadata_orig.rds"
ln -s $frangieh_proc_dir"grna/grna_assignments_ungrouped_metadata_control.rds" $sceptre2_data_dir"frangieh/control/grna/metadata_orig.rds"
ln -s $frangieh_proc_dir"grna/grna_assignments_ungrouped_metadata_ifn-gamma.rds" $sceptre2_data_dir"frangieh/ifn_gamma/grna/metadata_orig.rds"
# copy the protein metadata
ln -s $frangieh_proc_dir"protein/protein_expression_metadata_co-culture.rds" $sceptre2_data_dir"frangieh/co_culture/protein/metadata_orig.rds"
ln -s $frangieh_proc_dir"protein/protein_expression_metadata_control.rds" $sceptre2_data_dir"frangieh/control/protein/metadata_orig.rds"
ln -s $frangieh_proc_dir"protein/protein_expression_metadata_ifn-gamma.rds" $sceptre2_data_dir"frangieh/ifn_gamma/protein/metadata_orig.rds"
