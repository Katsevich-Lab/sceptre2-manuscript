# creates sym links to backing .odm files

source ~/.research_config
sceptre2_data_dir=$LOCAL_SCEPTRE2_DATA_DIR"data/"
papalexi_proc_dir=$LOCAL_PAPALEXI_2021_DATA_DIR"processed/"
liscovitch_proc_dir=$LOCAL_LISCOVITCH_2021_DATA_DIR"processed/"
schraivogel_proc_dir=$LOCAL_SCHRAIVOGEL_2020_DATA_DIR"processed/"
frangieh_proc_dir=$LOCAL_FRANGIEH_2021_DATA_DIR"processed/"

# Papalexi data
# gene
ln -sf $papalexi_proc_dir"gene/expression_matrix.odm" $sceptre2_data_dir"papalexi/eccite_screen/gene/matrix.odm"
ln -sf $papalexi_proc_dir"gene/metadata.rds" $sceptre2_data_dir"papalexi/eccite_screen/gene/metadata_orig.rds"
# grna_expression
ln -sf $papalexi_proc_dir"grna_expression/count_matrix.odm" $sceptre2_data_dir"papalexi/eccite_screen/grna_expression/matrix.odm"
ln -sf $papalexi_proc_dir"grna_expression/metadata.rds" $sceptre2_data_dir"papalexi/eccite_screen/grna_expression/metadata_orig.rds"
# grna_assignment
ln -sf $papalexi_proc_dir"grna_assignment/assignment_matrix.odm" $sceptre2_data_dir"papalexi/eccite_screen/grna_assignment/matrix.odm"
ln -sf $papalexi_proc_dir"grna_assignment/metadata.rds" $sceptre2_data_dir"papalexi/eccite_screen/grna_assignment/metadata_orig.rds"
# protein
ln -sf $papalexi_proc_dir"protein/count_matrix.odm" $sceptre2_data_dir"papalexi/eccite_screen/protein/matrix.odm"
ln -sf $papalexi_proc_dir"protein/metadata.rds" $sceptre2_data_dir"papalexi/eccite_screen/protein/metadata_orig.rds"

# Liscovitch data
# experiment_small chromatin
ln -sf $liscovitch_proc_dir"experiment_small/chromatin/chip_counts.odm" $sceptre2_data_dir"liscovitch/experiment_small/chromatin/matrix.odm"
ln -sf $liscovitch_proc_dir"experiment_small/chromatin/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_small/chromatin/metadata_orig.rds"
# experiment_small grna_expression
ln -sf $liscovitch_proc_dir"experiment_small/grna_expression/grna_counts.odm" $sceptre2_data_dir"liscovitch/experiment_small/grna_expression/matrix.odm"
ln -sf $liscovitch_proc_dir"experiment_small/grna_expression/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_small/grna_expression/metadata_orig.rds"
# experiment_small grna_assignment
ln -sf $liscovitch_proc_dir"experiment_small/grna_assignment/grna_assignments.odm" $sceptre2_data_dir"liscovitch/experiment_small/grna_assignment/matrix.odm"
ln -sf $liscovitch_proc_dir"experiment_small/grna_assignment/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_small/grna_assignment/metadata_orig.rds"
# experiment_big chromatin
ln -sf $liscovitch_proc_dir"experiment_big/chromatin/chip_counts.odm" $sceptre2_data_dir"liscovitch/experiment_big/chromatin/matrix.odm"
ln -sf $liscovitch_proc_dir"experiment_big/chromatin/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_big/chromatin/metadata_orig.rds"
# experiment_big grna_expression
ln -sf $liscovitch_proc_dir"experiment_big/grna_expression/grna_counts.odm" $sceptre2_data_dir"liscovitch/experiment_big/grna_expression/matrix.odm"
ln -sf $liscovitch_proc_dir"experiment_big/grna_expression/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_big/grna_expression/metadata_orig.rds"
# experiment_big grna_assignment
ln -sf $liscovitch_proc_dir"experiment_big/grna_assignment/grna_assignments.odm" $sceptre2_data_dir"liscovitch/experiment_big/grna_assignment/matrix.odm"
ln -sf $liscovitch_proc_dir"experiment_big/grna_assignment/metadata.rds" $sceptre2_data_dir"liscovitch/experiment_big/grna_assignment/metadata_orig.rds"

# Schraivogel data
# enhancer_screen_chr11 gene
ln -sf $schraivogel_proc_dir"enhancer_screen_chr11/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/gene/matrix.odm"
ln -sf $schraivogel_proc_dir"enhancer_screen_chr11/gene/metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/gene/metadata_orig.rds"
# enhancer_screen_chr11 grna_expression
ln -sf $schraivogel_proc_dir"enhancer_screen_chr11/grna_expression/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/grna_expression/matrix.odm"
ln -sf $schraivogel_proc_dir"enhancer_screen_chr11/grna_expression/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/grna_expression/metadata_orig.rds"
# enhancer_screen_chr11 grna_assignment
ln -sf $schraivogel_proc_dir"enhancer_screen_chr11/grna_assignment/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/grna_assignment/matrix.odm"
ln -sf $schraivogel_proc_dir"enhancer_screen_chr11/grna_assignment/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr11/grna_assignment/metadata_orig.rds"
# enhancer_screen_chr8 gene
ln -sf $schraivogel_proc_dir"enhancer_screen_chr8/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/gene/matrix.odm"
ln -sf $schraivogel_proc_dir"enhancer_screen_chr8/gene/metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/gene/metadata_orig.rds"
# enhancer_screen_chr8 grna_expression
ln -sf $schraivogel_proc_dir"enhancer_screen_chr8/grna_expression/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/grna_expression/matrix.odm"
ln -sf $schraivogel_proc_dir"enhancer_screen_chr8/grna_expression/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/grna_expression/metadata_orig.rds"
# enhancer_screen_chr8 grna_assignment
ln -sf $schraivogel_proc_dir"enhancer_screen_chr8/grna_assignment/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/grna_assignment/matrix.odm"
ln -sf $schraivogel_proc_dir"enhancer_screen_chr8/grna_assignment/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/enhancer_screen_chr8/grna_assignment/metadata_orig.rds"
# ground_truth_perturbseq gene
ln -sf $schraivogel_proc_dir"ground_truth_perturbseq/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/gene/matrix.odm"
ln -sf $schraivogel_proc_dir"ground_truth_perturbseq/gene/metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/gene/metadata_orig.rds"
# ground_truth_perturbseq grna_expression
ln -sf $schraivogel_proc_dir"ground_truth_perturbseq/grna_expression/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/grna_expression/matrix.odm"
ln -sf $schraivogel_proc_dir"ground_truth_perturbseq/grna_expression/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/grna_expression/metadata_orig.rds"
# ground_truth_perturbseq grna_assignment
ln -sf $schraivogel_proc_dir"ground_truth_perturbseq/grna_assignment/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/grna_assignment/matrix.odm"
ln -sf $schraivogel_proc_dir"ground_truth_perturbseq/grna_assignment/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_perturbseq/grna_assignment/metadata_orig.rds"
# ground_truth_tapseq gene
ln -sf $schraivogel_proc_dir"ground_truth_tapseq/gene/expression_matrix.odm" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/gene/matrix.odm"
ln -sf $schraivogel_proc_dir"ground_truth_tapseq/gene/metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/gene/metadata_orig.rds"
# ground_truth_tapseq grna_expression
ln -sf $schraivogel_proc_dir"ground_truth_tapseq/grna_expression/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/grna_expression/matrix.odm"
ln -sf $schraivogel_proc_dir"ground_truth_tapseq/grna_expression/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/grna_expression/metadata_orig.rds"
# ground_truth_tapseq grna_assignment
ln -sf $schraivogel_proc_dir"ground_truth_tapseq/grna_assignment/raw_ungrouped.odm" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/grna_assignment/matrix.odm"
ln -sf $schraivogel_proc_dir"ground_truth_tapseq/grna_assignment/raw_ungrouped_metadata.rds" $sceptre2_data_dir"schraivogel/ground_truth_tapseq/grna_assignment/metadata_orig.rds"

# Frangieh
# copy the gene odms
ln -sf $frangieh_proc_dir"co_culture/gene/gene_expression_matrix.odm" $sceptre2_data_dir"frangieh/co_culture/gene/matrix.odm"
ln -sf $frangieh_proc_dir"control/gene/gene_expression_matrix.odm" $sceptre2_data_dir"frangieh/control/gene/matrix.odm"
ln -sf $frangieh_proc_dir"ifn_gamma/gene/gene_expression_matrix.odm" $sceptre2_data_dir"frangieh/ifn_gamma/gene/matrix.odm"
# copy the grna_assignment odms
ln -sf $frangieh_proc_dir"co_culture/grna_assignment/grna_assignments_ungrouped.odm" $sceptre2_data_dir"frangieh/co_culture/grna_assignment/matrix.odm"
ln -sf $frangieh_proc_dir"control/grna_assignment/grna_assignments_ungrouped.odm" $sceptre2_data_dir"frangieh/control/grna_assignment/matrix.odm"
ln -sf $frangieh_proc_dir"ifn_gamma/grna_assignment/grna_assignments_ungrouped.odm" $sceptre2_data_dir"frangieh/ifn_gamma/grna_assignment/matrix.odm"
# copy the protein odms
ln -sf $frangieh_proc_dir"co_culture/protein/protein_expression_matrix.odm" $sceptre2_data_dir"frangieh/co_culture/protein/matrix.odm"
ln -sf $frangieh_proc_dir"control/protein/protein_expression_matrix.odm" $sceptre2_data_dir"frangieh/control/protein/matrix.odm"
ln -sf $frangieh_proc_dir"ifn_gamma/protein/protein_expression_matrix.odm" $sceptre2_data_dir"frangieh/ifn_gamma/protein/matrix.odm"
# copy the gene metadata
ln -sf $frangieh_proc_dir"co_culture/gene/gene_expression_metadata.rds" $sceptre2_data_dir"frangieh/co_culture/gene/metadata_orig.rds"
ln -sf $frangieh_proc_dir"control/gene/gene_expression_metadata.rds" $sceptre2_data_dir"frangieh/control/gene/metadata_orig.rds"
ln -sf $frangieh_proc_dir"ifn_gamma/gene/gene_expression_metadata.rds" $sceptre2_data_dir"frangieh/ifn_gamma/gene/metadata_orig.rds"
# copy the grna_assignment metadata
ln -sf $frangieh_proc_dir"co_culture/grna_assignment/grna_assignments_ungrouped_metadata.rds" $sceptre2_data_dir"frangieh/co_culture/grna_assignment/metadata_orig.rds"
ln -sf $frangieh_proc_dir"control/grna_assignment/grna_assignments_ungrouped_metadata.rds" $sceptre2_data_dir"frangieh/control/grna_assignment/metadata_orig.rds"
ln -sf $frangieh_proc_dir"ifn_gamma/grna_assignment/grna_assignments_ungrouped_metadata.rds" $sceptre2_data_dir"frangieh/ifn_gamma/grna_assignment/metadata_orig.rds"
# copy the protein metadata
ln -sf $frangieh_proc_dir"co_culture/protein/protein_expression_metadata.rds" $sceptre2_data_dir"frangieh/co_culture/protein/metadata_orig.rds"
ln -sf $frangieh_proc_dir"control/protein/protein_expression_metadata.rds" $sceptre2_data_dir"frangieh/control/protein/metadata_orig.rds"
ln -sf $frangieh_proc_dir"ifn_gamma/protein/protein_expression_metadata.rds" $sceptre2_data_dir"frangieh/ifn_gamma/protein/metadata_orig.rds"
