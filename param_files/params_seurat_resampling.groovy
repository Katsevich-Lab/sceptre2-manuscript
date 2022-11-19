data_method_pairs = ["frangieh/control/gene": ["mann_whitney_perm"]]

// FIRST, define the row and column names of the below matrices and vectors
row_names = ["frangieh/control/gene"]
col_names = ["mann_whitney_perm"]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
[5] // frangieh/control/gene
]
// mann_whitney_perm

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["all.q"]    // frangieh/control/gene
]
// mann_whitney_perm

// FOURTH, define an ordered list of optional arguments to each of the methods (Should be strings of the form "arg1=value1;arg2=value2;arg3=value3")
optional_args = [
"full_output=FALSE:B=200000" // mann_whitney_perm
]
