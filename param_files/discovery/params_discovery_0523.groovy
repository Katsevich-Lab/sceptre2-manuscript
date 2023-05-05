data_method_pairs_indiv = ["frangieh/control/gene": ["sceptre", "seurat_de", "liscovitch_method", "weissman_method", "seurat_de_nb", "mimosca", "schraivogel_method"],
                           "papalexi/eccite_screen/gene": ["sceptre", "seurat_de", "liscovitch_method", "weissman_method", "seurat_de_nb", "mimosca", "schraivogel_method"]
                          ]
                            
data_method_pairs_grouped = []

row_names = ["frangieh/control/gene",
             "papalexi/eccite_screen/gene"]
             
col_names = ["schraivogel_method",
             "seurat_de",
             "liscovitch_method",
             "mimosca",
             "weissman_method",
             "sceptre",
             "seurat_de_nb"]
             
// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
                         [34, 13, 13, 35, 49, 8, 13], // frangieh/co_culture/gene
                         [23,  9, 45, 26, 34, 8, 9], // frangieh/control/gene
                         [36, 13, 13, 34, 53, 8, 13], // frangieh/ifn_gamma/gene
                         [18,  6,  6, 19, 22, 6, 6], // papalexi/eccite_screen/gene
                         [2,   1,  1,  5,  1, 4, 1], // papalexi/eccite_screen/protein
                         [6,  12, 12, 18,  6, 8, 12], // schraivogel/enhancer_screen_chr11/gene
                         [8,  14, 14, 21,  6, 8, 14], // schraivogel/enhancer_screen_chr8/gene
                         ]
// schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method, sceptre, seurat_de_nb

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
["all.q", "all.q", "all.q", "all.q", "all.q", "all.q", "all.q"], // frangieh/control/gene
["all.q", "all.q", "all.q", "all.q", "all.q", "all.q", "all.q"]  // papalexi/eccite_screen/gene
]
// schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method, sceptre, seurat_de_nb

// FOURTH, define an ordered list of optional arguments to each of the methods (Should be strings of the form "arg1=value1;arg2=value2;arg3=value3")
optional_args = [
"", // schraivogel_method
"", // seurat_de
"", // liscovitch_method
"", // mimosca
"", // weissman_method
"", // sceptre
""  // seurat_de_nb
]