// FIRST, define the dataset-method pairs to analyze in a map
/*
data_method_pairs = ["frangieh/co_culture/gene": ["schraivogel_method", "seurat_de", "mimosca", "weissman_method"],
                     "frangieh/control/gene": ["schraivogel_method", "seurat_de", "mimosca", "weissman_method"],
                     "frangieh/ifn_gamma/gene": ["schraivogel_method", "seurat_de", "mimosca", "weissman_method"],
                     "liscovitch/experiment_big/chromatin": ["liscovitch_method"],
                     "liscovitch/experiment_small/chromatin": ["liscovitch_method"],
                     "papalexi/eccite_screen/gene": ["schraivogel_method", "seurat_de", "mimosca", "weissman_method"],
                     "papalexi/eccite_screen/protein": ["schraivogel_method", "seurat_de", "mimosca", "weissman_method"],
                     "schraivogel/enhancer_screen_chr11/gene": ["schraivogel_method", "seurat_de", "mimosca", "weissman_method"],
                     "schraivogel/enhancer_screen_chr8/gene": ["schraivogel_method", "seurat_de", "mimosca", "weissman_method"],
                     "schraivogel/ground_truth_perturbseq/gene": ["schraivogel_method", "seurat_de", "mimosca", "weissman_method"],
                     "schraivogel/ground_truth_tapseq/gene": ["schraivogel_method", "seurat_de", "mimosca", "weissman_method"],
                     "simulated/experiment_1/gene": ["schraivogel_method", "seurat_de", "mimosca", "weissman_method"]
                     ]
*/
data_method_pairs = ["frangieh/co_culture/gene": ["schraivogel_method", "weissman_method"],
                     "frangieh/control/gene": ["schraivogel_method", "weissman_method"],
                     "frangieh/ifn_gamma/gene": ["schraivogel_method", "weissman_method"],
                     "papalexi/eccite_screen/gene": ["schraivogel_method", "weissman_method"],
                     "papalexi/eccite_screen/protein": ["schraivogel_method", "weissman_method"],
                     "schraivogel/enhancer_screen_chr11/gene": ["schraivogel_method", "weissman_method"],
                     "schraivogel/enhancer_screen_chr8/gene": ["schraivogel_method", "weissman_method"],
                     "schraivogel/ground_truth_perturbseq/gene": ["schraivogel_method", "weissman_method"],
                     "schraivogel/ground_truth_tapseq/gene": ["schraivogel_method", "weissman_method"],
                     "simulated/experiment_1/gene": ["schraivogel_method", "weissman_method"]
                     ]


// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
[34, 13, 1, 31, 45],
[23, 9, 1, 22, 45],
[36, 13, 1, 29, 45],
[2, 1, 1, 1, 45],
[2, 1, 1, 1, 45],
[18, 6, 1, 14, 45],
[2, 1, 1, 1, 45],
[6, 12, 1, 13, 45],
[8, 14, 1, 16, 45],
[44, 11, 1, 26, 45],
[2, 1, 1, 1, 45],
[28, 7, 1, 12, 45]
]


// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["short.q", "short.q", "short.q", "all.q", "short.q"],
["short.q", "short.q", "short.q", "all.q", "short.q"],
["short.q", "short.q", "short.q", "all.q", "short.q"],
["short.q", "short.q", "short.q", "short.q", "short.q"],
["short.q", "short.q", "short.q", "short.q", "short.q"],
["short.q", "short.q", "short.q", "all.q", "short.q"],
["short.q", "short.q", "short.q", "short.q", "short.q"],
["short.q", "short.q", "short.q", "short.q", "short.q"],
["short.q", "short.q", "short.q", "short.q", "short.q"],
["short.q", "short.q", "short.q", "all.q", "short.q"],
["short.q", "short.q", "short.q", "short.q", "short.q"],
["short.q", "short.q", "short.q", "all.q", "short.q"]
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
             "mimosca",
             "weissman_method"]


// FIFTH, define an ordered list of optional arguments to each of the methods
// Should be strings of the form "arg1=value1;arg2=value2;arg3=value3".
// If no optional argument, then set NA
optional_args = [
"", // schraivogel_method
"", // seurat_de
"", // liscovitch_method
"n_rep=10", // mimosca
"" // weissman_method
]
