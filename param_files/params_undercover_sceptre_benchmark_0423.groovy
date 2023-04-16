data_method_pairs = ["frangieh/co_culture/gene": ["sceptre_approximate_poisson", "sceptre_exact_poisson", "sceptre_approximate_nb", "sceptre_exact_nb"],
                     "frangieh/control/gene": ["sceptre_approximate_poisson", "sceptre_exact_poisson", "sceptre_approximate_nb", "sceptre_exact_nb"],
                     "frangieh/ifn_gamma/gene": ["sceptre_approximate_poisson", "sceptre_exact_poisson", "sceptre_approximate_nb", "sceptre_exact_nb"],
                     "papalexi/eccite_screen/gene": ["sceptre_approximate_poisson", "sceptre_exact_poisson", "sceptre_approximate_nb", "sceptre_exact_nb"],
                     "papalexi/eccite_screen/protein": ["sceptre_approximate_poisson", "sceptre_exact_poisson", "sceptre_approximate_nb", "sceptre_exact_nb"],
                     "schraivogel/enhancer_screen_chr11/gene": ["sceptre_approximate_poisson", "sceptre_exact_poisson", "sceptre_approximate_nb", "sceptre_exact_nb"],
                     "schraivogel/enhancer_screen_chr8/gene": ["sceptre_approximate_poisson", "sceptre_exact_poisson", "sceptre_approximate_nb", "sceptre_exact_nb"],
                     "simulated/experiment_1/gene": ["sceptre_approximate_poisson", "sceptre_exact_poisson", "sceptre_approximate_nb", "sceptre_exact_nb"]
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
col_names = ["sceptre_approximate_poisson",
             "sceptre_exact_poisson",
             "sceptre_approximate_nb",
             "sceptre_exact_nb"]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
[x, 8, 8, 8], // frangieh/co_culture/gene
[8, 8, 8, 8], // frangieh/control/gene
[8, 8, 8, 8], // frangieh/ifn_gamma/gene
[6, 6, 6, 6], // papalexi/eccite_screen/gene
[4, 4, 4, 4], // papalexi/eccite_screen/protein
[8, 8, 8, 8], // schraivogel/enhancer_screen_chr11/gene
[8, 8, 8, 8], // schraivogel/enhancer_screen_chr8/gene
[10, 10, 10, 10]  // simulated/experiment_1/gene
]
// sceptre_approximate_poisson, sceptre_exact_poisson, sceptre_approximate_nb, sceptre_exact_nb

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["short.q", "short.q", "short.q", "short.q"], // frangieh/co_culture/gene
["short.q", "short.q", "short.q", "short.q"], // frangieh/control/gene
["short.q", "short.q", "short.q", "short.q"], // frangieh/ifn_gamma/gene
["short.q", "short.q", "short.q", "short.q"], // papalexi/eccite_screen/gene
["short.q", "short.q", "short.q", "short.q"], // papalexi/eccite_screen/protein
["short.q", "short.q", "short.q", "short.q"], // schraivogel/enhancer_screen_chr11/gene
["short.q", "short.q", "short.q", "short.q"], // schraivogel/enhancer_screen_chr8/gene
["short.q", "short.q", "short.q", "short.q"]  // simulated/experiment_1/gene
]
// sceptre_approximate_poisson, sceptre_exact_poisson, sceptre_approximate_nb, sceptre_exact_nb

// FOURTH, define an ordered list of optional arguments to each of the methods (Should be strings of the form "arg1=value1;arg2=value2;arg3=value3")
optional_args = [
"", // sceptre_approximate_poisson
"", // sceptre_exact_poisson
"", // sceptre_approximate_nb
"", // sceptre_exact_nb
]
