// FIRST, define the dataset-method pairs to analyze in a map

/*
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
*/
data_method_pairs = ["frangieh/co_culture/gene": ["mimosca"],
                     "frangieh/control/gene": ["mimosca"],
                     "frangieh/ifn_gamma/gene": ["mimosca"],
                     "papalexi/eccite_screen/gene": ["mimosca"],
                     "papalexi/eccite_screen/protein": ["mimosca"],
                     "schraivogel/enhancer_screen_chr11/gene": ["mimosca"],
                     "schraivogel/enhancer_screen_chr8/gene": ["mimosca"],
                     "schraivogel/ground_truth_perturbseq/gene": ["mimosca"],
                     "schraivogel/ground_truth_tapseq/gene": ["mimosca"],
                     "simulated/experiment_1/gene": ["mimosca"]]


// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
[34, 13, 1, 40],
[23, 9, 1, 40],
[36, 13, 1, 40],
[2, 1, 1, 40],
[2, 1, 1, 40],
[18, 6, 1, 40],
[2, 1, 1, 40],
[6, 12, 1, 40],
[8, 14, 1, 40],
[44, 11, 1, 40],
[2, 1, 1, 40],
[28, 7, 1, 40]
]


// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["-q short.q", "-q short.q", "-q short.q", "-q all.q"],
["-q short.q", "-q short.q", "-q short.q", "-q all.q"],
["-q short.q", "-q short.q", "-q short.q", "-q all.q"],
["-q short.q", "-q short.q", "-q short.q", "-q short.q"],
["-q short.q", "-q short.q", "-q short.q", "-q short.q"],
["-q short.q", "-q short.q", "-q short.q", "-q all.q"],
["-q short.q", "-q short.q", "-q short.q", "-q short.q"],
["-q short.q", "-q short.q", "-q short.q", "-q short.q"],
["-q short.q", "-q short.q", "-q short.q", "-q short.q"],
["-q short.q", "-q short.q", "-q short.q", "-q all.q"],
["-q short.q", "-q short.q", "-q short.q", "-q short.q"],
["-q short.q", "-q short.q", "-q short.q", "-q all.q"]
]


// FOURTH, define the row and column names of the above matrices
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
             "seurat_de",
             "liscovitch_method",
             "mimosca"]
