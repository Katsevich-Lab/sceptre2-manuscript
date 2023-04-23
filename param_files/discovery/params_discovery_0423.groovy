data_method_pairs_indiv = ["frangieh/control/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method"],
                           "papalexi/eccite_screen/gene": ["schraivogel_method", "seurat_de", "liscovitch_method", "mimosca", "weissman_method"]
                          ]
                            
data_method_pairs_grouped = []

row_names = ["frangieh/control/gene",
             "papalexi/eccite_screen/gene"]
             
col_names = ["sceptre",
             "schraivogel_method",
             "seurat_de",
             "liscovitch_method",
             "mimosca",
             "weissman_method"]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
[8, 23,  9, 45, 26, 34], // frangieh/control/gene
[6, 18,  6,  6, 19, 22]  // papalexi/eccite_screen/gene
]
// sceptre, schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["all.q", "all.q", "all.q", "all.q", "all.q", "all.q"], // frangieh/control/gene
["all.q", "all.q", "all.q", "all.q", "all.q", "all.q"]  // papalexi/eccite_screen/gene
]
// sceptre, schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method

// FOURTH, define an ordered list of optional arguments to each of the methods (Should be strings of the form "arg1=value1;arg2=value2;arg3=value3")
optional_args = [
"", // sceptre
"", // schraivogel_method
"", // seurat_de
"", // liscovitch_method
"", // mimosca
"", // weissman_method
]