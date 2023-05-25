data_method_pairs = ["frangieh/ifn_gamma/gene": ["sceptre_orig"],
                     "papalexi/eccite_screen/gene": ["weissman_method"],
                     "papalexi/eccite_screen/protein": ["weissman_method"]
                     ]
                     
row_names = ["frangieh/ifn_gamma/gene",
             "papalexi/eccite_screen/gene",
             "papalexi/eccite_screen/protein"]
col_names = ["sceptre_orig",
             "weissman_method"]

data_method_ram_matrix = [
[30, 53], // frangieh/ifn_gamma/gene
[15, 22],  // papalexi/eccite_screen/gene
[15, 5]   // papalexi/eccite_screen/protein
]
// sceptre_orig, weissman_method

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["all.q", "short.q"], // frangieh/ifn_gamma/gene
["all.q", "short.q"], // papalexi/eccite_screen/gene
["all.q", "short.q"], // papalexi/eccite_screen/protein
]
// sceptre_orig, weissman_method

optional_args = [
"", // sceptre_orig
"use_batch=FALSE"  // weissman_method
]