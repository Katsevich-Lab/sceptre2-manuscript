// FIRST, define the dataset-method pairs to analyze in a map
data_method_pairs = ["frangieh/co_culture/gene": ["schraivogel_method", "seurat_de"],
                     "frangieh/control/gene": ["schraivogel_method", "seurat_de"],
                     "frangieh/ifn_gamma/gene": ["schraivogel_method", "seurat_de"],
                     "liscovitch/experiment_big/chromatin": ["schraivogel_method", "seurat_de"],
                     "liscovitch/experiment_small/chromatin": ["schraivogel_method", "seurat_de"],
                     "papalexi/eccite_screen/gene": ["schraivogel_method", "seurat_de"],
                     "papalexi/eccite_screen/protein": ["schraivogel_method", "seurat_de"],
                     "schraivogel/enhancer_screen_chr11/gene": ["schraivogel_method", "seurat_de"],
                     "schraivogel/enhancer_screen_chr8/gene": ["schraivogel_method", "seurat_de"],
                     "schraivogel/ground_truth_perturbseq/gene": ["schraivogel_method", "seurat_de"],
                     "schraivogel/ground_truth_tapseq/gene": ["schraivogel_method", "seurat_de"],
                     "simulated/experiment_1/gene": ["schraivogel_method", "seurat_de"]]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [[34, 13],
[23, 9],
[36, 13],
[2, 1],
[2, 1],
[18, 6],
[2, 1],
[6, 12],
[8, 14],
[44, 11],
[2, 1],
[28, 7]]

row_names = ["frangieh/co_culture/gene",
             "frangieh/control/gene",
             "frangieh/ifn_gamma/gene",
             "liscovitch/experiment_big/chromatin",
             "liscovitch/experiment_small/chromatin",
             "papalexi/eccite_screen/gene",
             "papalexi/eccite_screen/protein",
             "schraivogel/enhancer_screen_chr11/gene",
             "schraivogel/enhancer_screen_chr8/gene",
             "schraivogel/ground_truth_perturbseq/gene",
             "schraivogel/ground_truth_tapseq/gene",
             "simulated/experiment_1/gene"]

col_names = ["schraivogel_method",
             "seurat_de"]
