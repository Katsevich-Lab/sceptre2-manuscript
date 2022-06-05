// FIRST, define the dataset-method pairs to analyze in a map
data_method_pairs = [
 "frangieh/co_culture/gene": ["dummy_method_1", "dummy_method_2"],
 "papalexi/eccite_screen/protein": ["dummy_method_2"],
 "schraivogel/ground_truth_tapseq/gene": ["dummy_method_1", "dummy_method_2"],
 "simulated/experiment_1/gene": ["dummy_method_1"]
]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
 [1, 1],
 [1, 1],
 [1, 1],
 [1, 1],
 [1, 1],
 [1, 1],
 [1, 1],
 [1, 1],
 [1, 1],
 [1, 1],
 [1, 1],
 [1, 1]
 ]


 // THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
 ["short.q", "short.q"],
 ["short.q", "short.q"],
 ["short.q", "short.q"],
 ["short.q", "short.q"],
 ["short.q", "short.q"],
 ["short.q", "short.q"],
 ["short.q", "short.q"],
 ["short.q", "short.q"],
 ["short.q", "short.q"],
 ["short.q", "short.q"],
 ["short.q", "short.q"],
 ["short.q", "short.q"]
]


// FOURTH, define the row and column names of the above matrices
row_names = [
 "frangieh/co_culture/gene",
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
col_names = [
"dummy_method_1",
"dummy_method_2"
]


// FIFTH, define an list of optional arguments (where the order of methods corresponds to col_names, above)
// Should be strings of the form "arg1=value1 arg2=value2 arg3=value3".
// If no optional argument, then set NA
optional_args = [
 "dist=rnorm",
 "shape1=10 shape2=1"
]
