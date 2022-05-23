// FIRST, define the dataset-method pairs to analyze in a map

data_method_pairs = ["schraivogel/ground_truth_tapseq/gene": ["schraivogel_method", "seurat_de"],
                     "schraivogel/ground_truth_perturbseq/gene": ["schraivogel_method", "seurat_de"],
                     "schraivogel/enhancer_screen_chr11/gene": ["schraivogel_method", "seurat_de"],
                     "schraivogel/enhancer_screen_chr8/gene": ["schraivogel_method", "seurat_de"],
                     "papalexi/eccite_screen/gene": ["schraivogel_method", "seurat_de"],
                     "papalexi/eccite_screen/protein": ["schraivogel_method", "seurat_de"],
                     "frangieh/co_culture/gene": ["schraivogel_method", "seurat_de"],
                     "frangieh/control/gene": ["schraivogel_method", "seurat_de"],
                     "frangieh/ifn_gamma/gene": ["schraivogel_method", "seurat_de"],
                     "liscovitch/experiment_small/chromatin": ["schraivogel_method", "seurat_de"],
                     "liscovitch/experiment_big/chromatin": ["schraivogel_method", "seurat_de"],
                     "simulated/experiment_1/gene": ["schraivogel_method", "seurat_de"]]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [[ 8,  8, 1, 1],
                          [45, 20, 1, 1],
                          [45, 20, 1, 1],
                          [45, 20, 1, 1],
                          [45, 20, 1, 1],
                          [45, 20, 1, 1],
                          [45, 20, 1, 1],
                          [45, 20, 1, 1],
                          [45, 20, 1, 1],
                          [45, 20, 1, 1],
                          [45, 20, 1, 1],
                          [45, 20, 1, 1]]

row_names = ["schraivogel/ground_truth_tapseq/gene",
             "schraivogel/ground_truth_perturbseq/gene",
             "schraivogel/enhancer_screen_chr11/gene",
             "schraivogel/enhancer_screen_chr8/gene",
             "papalexi/eccite_screen/gene",
             "papalexi/eccite_screen/protein",
             "frangieh/co_culture/gene",
             "frangieh/control/gene",
             "frangieh/ifn_gamma/gene",
             "liscovitch/experiment_small/chromatin",
             "liscovitch/experiment_big/chromatin",
             "simulated/experiment_1/gene"]

col_names = ["schraivogel_method",
             "seurat_de",
             "dummy_method_1",
             "dummy_method_2"]
