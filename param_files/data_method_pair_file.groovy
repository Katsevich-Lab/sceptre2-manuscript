// FIRST, define the dataset-method pairs to analyze in a map
/*
data_method_pairs = ["frangieh/co_culture/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method", "nb_regression", "fisher_exact", "sceptre"],
                     "frangieh/control/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method", "nb_regression", "fisher_exact", "sceptre"],
                     "frangieh/ifn_gamma/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method", "nb_regression", "fisher_exact", "sceptre"],
                     "liscovitch/experiment_big/chromatin": ["liscovitch_method"],
                     "liscovitch/experiment_small/chromatin": ["liscovitch_method"],
                     "papalexi/eccite_screen/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method", "nb_regression", "fisher_exact", "sceptre"],
                     "papalexi/eccite_screen/protein": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method", "nb_regression", "fisher_exact", "sceptre"],
                     "schraivogel/enhancer_screen_chr11/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method", "nb_regression", "fisher_exact", "sceptre"],
                     "schraivogel/enhancer_screen_chr8/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method", "nb_regression", "fisher_exact", "sceptre"],
                     "schraivogel/ground_truth_perturbseq/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method", "nb_regression", "fisher_exact", "sceptre"],
                     "schraivogel/ground_truth_tapseq/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method", "nb_regression", "fisher_exact", "sceptre"],
                     "simulated/experiment_1/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method", "nb_regression", "fisher_exact", "sceptre"]
                     ]
*/

data_method_pairs = ["schraivogel/ground_truth_perturbseq/gene": ["fisher_exact"],
                     "schraivogel/ground_truth_tapseq/gene": ["fisher_exact"]
                     ]

// FIRST, define the row and column names of the below matrices and vectors
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
             "weissman_method",
             "nb_regression",
             "fisher_exact",
             "sceptre"]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
[34, 13, 13, 35, 49, 5, 5, 5], // frangieh/co_culture/gene
[23,  9, 45, 26, 34, 5, 5, 5], // frangieh/control/gene
[36, 13, 13, 34, 53, 5, 5, 5], // frangieh/ifn_gamma/gene
[2,   1,  1,  1,  1, 1, 1, 1], // liscovitch/experiment_big/chromatin
[2,   1,  1,  1,  1, 1, 1, 1], // liscovitch/experiment_small/chromatin
[18,  6,  6, 19, 22, 5, 5, 5], // papalexi/eccite_screen/gene
[2,   1,  1,  5,  1, 2, 5, 5], // papalexi/eccite_screen/protein
[6,  12, 12, 18,  6, 5, 5, 5], // schraivogel/enhancer_screen_chr11/gene
[8,  14, 14, 21,  6, 5, 5, 5], // schraivogel/enhancer_screen_chr8/gene
[44, 11, 11, 31, 34, 5, 5, 5], // schraivogel/ground_truth_perturbseq/gene
[2,   1,  1,  5,  1, 2, 5, 5], // schraivogel/ground_truth_tapseq/gene
[28,  7,  7, 17, 19, 5, 5, 5] // simulated/experiment_1/gene
]
// schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method, nb_regression, fisher_exact, sceptre

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["short.q", "short.q", "short.q", "all.q",   "short.q", "all.q",   "short.q", "all.q"],  // frangieh/co_culture/gene
["short.q", "short.q", "short.q", "all.q",   "short.q", "all.q",   "short.q", "all.q"], // frangieh/control/gene
["short.q", "short.q", "short.q", "all.q",   "short.q", "all.q",   "short.q", "all.q"],  // frangieh/ifn_gamma/gene
["short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q"], // liscovitch/experiment_big/chromatin
["short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q"], // liscovitch/experiment_small/chromatin
["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q", "all.q"],  // papalexi/eccite_screen/gene
["short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q"], // papalexi/eccite_screen/protein
["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q", "short.q"], // schraivogel/enhancer_screen_chr11/gene
["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q", "short.q"],  // schraivogel/enhancer_screen_chr8/gene
["short.q", "short.q", "short.q", "all.q",   "short.q", "all.q",   "short.q", "all.q"], // schraivogel/ground_truth_perturbseq/gene
["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q", "short.q"], // schraivogel/ground_truth_tapseq/gene
["short.q", "short.q", "short.q", "all.q",   "short.q",   "all.q",   "short.q", "all.q"] // simulated/experiment_1/gene
]
// schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method, nb_regression, fisher_exact, sceptre

// FOURTH, define an ordered list of optional arguments to each of the methods (Should be strings of the form "arg1=value1;arg2=value2;arg3=value3")
optional_args = [
"", // schraivogel_method
"", // seurat_de
"", // liscovitch_method
"n_rep=50", // mimosca
"", // weissman_method
"progress=FALSE", // nb_regression
"progress=FALSE;control_group=compliment", // fisher_exact
"output_amount=1" // sceptre
]
