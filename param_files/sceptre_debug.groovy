// FIRST, define the dataset-method pairs to analyze in a map
data_method_pairs = ["frangieh/ifn_gamma/gene": ["sceptre"],
                     "papalexi/eccite_screen/gene": ["sceptre"],
                     "papalexi/eccite_screen/protein": ["sceptre"]
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
col_names = ["sceptre"]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
[10], // frangieh/co_culture/gene
[10], // frangieh/control/gene
[10], // frangieh/ifn_gamma/gene
[1], // liscovitch/experiment_big/chromatin
[1], // liscovitch/experiment_small/chromatin
[10], // papalexi/eccite_screen/gene
[10], // papalexi/eccite_screen/protein
[5], // schraivogel/enhancer_screen_chr11/gene
[5], // schraivogel/enhancer_screen_chr8/gene
[5], // schraivogel/ground_truth_perturbseq/gene
[5], // schraivogel/ground_truth_tapseq/gene
[5] // simulated/experiment_1/gene
]
// schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method, permutation_test, nb_regression, sceptre

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["all.q"],  // frangieh/co_culture/gene
["all.q"], // frangieh/control/gene
["all.q"],  // frangieh/ifn_gamma/gene
["short.q"], // liscovitch/experiment_big/chromatin
["short.q"], // liscovitch/experiment_small/chromatin
["all.q"],  // papalexi/eccite_screen/gene
["short.q"], // papalexi/eccite_screen/protein
["short.q"], // schraivogel/enhancer_screen_chr11/gene
["short.q"],  // schraivogel/enhancer_screen_chr8/gene
["all.q"], // schraivogel/ground_truth_perturbseq/gene
["short.q"], // schraivogel/ground_truth_tapseq/gene
["all.q"] // simulated/experiment_1/gene
]
// schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method, permutation_test, nb_regression, sceptre

// FOURTH, define an ordered list of optional arguments to each of the methods (Should be strings of the form "arg1=value1:arg2=value2:arg3=value3")
optional_args = [
"output_amount=1:B=300000:sn_approx=FALSE" // sceptre
]
