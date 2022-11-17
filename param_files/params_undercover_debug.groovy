data_method_pairs = ["frangieh/co_culture/gene": ["nb_regression_w_covariates", "nb_regression_no_covariates", "sceptre"],
                     "frangieh/control/gene": ["nb_regression_w_covariates", "nb_regression_no_covariates", "sceptre"],
                     "frangieh/ifn_gamma/gene": ["nb_regression_w_covariates", "nb_regression_no_covariates", "sceptre"],
                     "papalexi/eccite_screen/gene": ["nb_regression_w_covariates", "nb_regression_no_covariates", "sceptre"],
                     "papalexi/eccite_screen/protein": ["nb_regression_w_covariates", "nb_regression_no_covariates", "sceptre"],
                     "schraivogel/enhancer_screen_chr11/gene": ["nb_regression_w_covariates", "nb_regression_no_covariates", "sceptre"],
                     "schraivogel/enhancer_screen_chr8/gene": ["nb_regression_w_covariates", "nb_regression_no_covariates", "sceptre"],
                     "simulated/experiment_1/gene": ["nb_regression_w_covariates", "nb_regression_no_covariates", "sceptre"]
                     ]

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
             "sceptre"]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
[5, 5, 5], // frangieh/co_culture/gene
[5, 5, 5], // frangieh/control/gene
[5, 5, 5], // frangieh/ifn_gamma/gene
[5, 5, 6], // papalexi/eccite_screen/gene
[5, 5, 4], // papalexi/eccite_screen/protein
[5, 5, 9], // schraivogel/enhancer_screen_chr11/gene
[5, 5, 9], // schraivogel/enhancer_screen_chr8/gene
[5, 5, 7]  // simulated/experiment_1/gene
]
// nb_regression_w_covariates, nb_regression_no_covariates, sceptre

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["all.q",   "all.q",   "all.q"],   // frangieh/co_culture/gene
["all.q",   "all.q",   "all.q"],   // frangieh/control/gene
["all.q",   "all.q",   "all.q"],   // frangieh/ifn_gamma/gene
["all.q",   "all.q",   "all.q"],   // papalexi/eccite_screen/gene
["short.q", "short.q", "short.q"], // papalexi/eccite_screen/protein
["short.q", "short.q", "short.q"],   // schraivogel/enhancer_screen_chr11/gene
["short.q", "short.q", "short.q"],   // schraivogel/enhancer_screen_chr8/gene
["short.q", "short.q", "short.q"]    // simulated/experiment_1/gene
]
// nb_regression_w_covariates, nb_regression_no_covariates, sceptre

// FOURTH, define an ordered list of optional arguments to each of the methods (Should be strings of the form "arg1=value1;arg2=value2;arg3=value3")
optional_args = [
"", // nb_regression_w_covariates
"", // nb_regression_no_covariates
"output_amount=1:B=2500000:with_covariates=FALSE" // sceptre
]
