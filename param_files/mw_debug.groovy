// FIRST, define the dataset-method pairs to analyze in a map
data_method_pairs = ["frangieh/co_culture/gene": ["mann_whitney_perm_mimosca"],
                     "frangieh/control/gene": ["mann_whitney_perm_mimosca"],
                     "frangieh/ifn_gamma/gene": ["mann_whitney_perm_mimosca"],
                     "papalexi/eccite_screen/gene": ["mann_whitney_perm_mimosca"],
                     "papalexi/eccite_screen/protein": ["mann_whitney_perm_mimosca"],
                     "schraivogel/ground_truth_perturbseq/gene": ["mann_whitney_perm_mimosca"],
                     "schraivogel/ground_truth_tapseq/gene": ["mann_whitney_perm_mimosca"],
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
col_names = ["mann_whitney_perm_mimosca"]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
[13], // frangieh/co_culture/gene
[45], // frangieh/control/gene
[13], // frangieh/ifn_gamma/gene
[1], // liscovitch/experiment_big/chromatin
[1], // liscovitch/experiment_small/chromatin
[6], // papalexi/eccite_screen/gene
[1], // papalexi/eccite_screen/protein
[12], // schraivogel/enhancer_screen_chr11/gene
[14], // schraivogel/enhancer_screen_chr8/gene
[11], // schraivogel/ground_truth_perturbseq/gene
[1], // schraivogel/ground_truth_tapseq/gene
[7] // simulated/experiment_1/gene
]
// mann_whitney_perm

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["short.q"],  // frangieh/co_culture/gene
["short.q"], // frangieh/control/gene
["short.q"],  // frangieh/ifn_gamma/gene
["short.q"], // liscovitch/experiment_big/chromatin
["short.q"], // liscovitch/experiment_small/chromatin
["short.q"],  // papalexi/eccite_screen/gene
["short.q"], // papalexi/eccite_screen/protein
["short.q"], // schraivogel/enhancer_screen_chr11/gene
["short.q"],  // schraivogel/enhancer_screen_chr8/gene
["short.q"], // schraivogel/ground_truth_perturbseq/gene
["short.q"], // schraivogel/ground_truth_tapseq/gene
["short.q"] // simulated/experiment_1/gene
]
// mann_whitney_perm

// FOURTH, define an ordered list of optional arguments to each of the methods (Should be strings of the form "arg1=value1:arg2=value2:arg3=value3")
optional_args = [
  "B=1" // mann_whitney_perm
]
