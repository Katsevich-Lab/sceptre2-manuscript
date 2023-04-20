data_method_pairs = ["papalexi/eccite_screen/gene": ["sceptre"]]
                     
// FIRST, define the row and column names of the below matrices and vectors
row_names = ["frangieh/co_culture/gene",
             "frangieh/control/gene",
             "frangieh/ifn_gamma/gene",
             "papalexi/eccite_screen/gene",
             "papalexi/eccite_screen/protein",
             "schraivogel/enhancer_screen_chr11/gene",
             "schraivogel/enhancer_screen_chr8/gene",
             "simulated/experiment_1/gene"]
col_names = ["nb_regression_w_covariates",
             "nb_regression_no_covariates",
             "sceptre",
             "sceptre_no_covariates",
             "schraivogel_method",
             "seurat_de",
             "liscovitch_method",
             "mimosca",
             "weissman_method"]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
[5, 5, 8, 8, 34, 13, 13, 35, 49], // frangieh/co_culture/gene
[5, 5, 8, 8, 23,  9, 45, 26, 34], // frangieh/control/gene
[5, 5, 8, 8, 36, 13, 13, 34, 53], // frangieh/ifn_gamma/gene
[5, 5, 6, 6, 18,  6,  6, 19, 22], // papalexi/eccite_screen/gene
[5, 5, 4, 4,  5,  5,  5,  5,  5], // papalexi/eccite_screen/protein
[5, 5, 8, 8,  6, 12, 12, 18,  6], // schraivogel/enhancer_screen_chr11/gene
[5, 5, 8, 8,  8, 14, 14, 21,  6], // schraivogel/enhancer_screen_chr8/gene
[5, 5, 10, 10, 28,  7,  7, 17, 19]  // simulated/experiment_1/gene
]
// nb_regression_w_covariates, nb_regression_no_covariates, sceptre, sceptre_no_covariates, schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["all.q",   "all.q",   "short.q", "short.q", "short.q", "short.q", "short.q", "all.q", "short.q"], // frangieh/co_culture/gene
["all.q",   "all.q",   "short.q", "short.q", "short.q", "short.q", "short.q", "all.q", "short.q"], // frangieh/control/gene
["all.q",   "all.q",   "short.q", "short.q", "short.q", "short.q", "short.q", "all.q", "short.q"], // frangieh/ifn_gamma/gene
["all.q",   "all.q",   "short.q", "short.q", "short.q", "short.q", "short.q", "all.q", "short.q"], // papalexi/eccite_screen/gene
["short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "all.q", "short.q"], // papalexi/eccite_screen/protein
["short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "all.q", "short.q"], // schraivogel/enhancer_screen_chr11/gene
["short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "all.q", "short.q"], // schraivogel/enhancer_screen_chr8/gene
["short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "all.q", "short.q"]  // simulated/experiment_1/gene
]
// nb_regression_w_covariates, nb_regression_no_covariates, sceptre, sceptre_no_covariates, schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method

// FOURTH, define an ordered list of optional arguments to each of the methods (Should be strings of the form "arg1=value1;arg2=value2;arg3=value3")
optional_args = [
"", // nb_regression_w_covariates
"", // nb_regression_no_covariates
"", // sceptre
"", // sceptre_no_covariates
"", // schraivogel_method
"", // seurat_de
"", // liscovitch_method
"", // mimosca
"", // weissman_method
]
