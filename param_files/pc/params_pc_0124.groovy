/*
data_method_pairs_indiv = ["frangieh/co_culture/gene": ["mimosca", "schraivogel_method"],
                           "frangieh/control/gene": ["mimosca", "schraivogel_method"],
                           "frangieh/ifn_gamma/gene": ["mimosca", "schraivogel_method"],
                           "papalexi/eccite_screen/gene": ["mimosca", "schraivogel_method"],
                           "papalexi/eccite_screen/protein": ["mimosca", "schraivogel_method"],
                           "schraivogel/enhancer_screen_chr11/gene": ["mimosca", "schraivogel_method"],
                           "schraivogel/enhancer_screen_chr8/gene": ["mimosca", "schraivogel_method"],
                           "simulated/experiment_2/gene": ["mimosca", "schraivogel_method"]
]
data_method_pairs_grouped = ["frangieh/co_culture/gene": ["seurat_de", "liscovitch_method", "weissman_method", "seurat_de_nb"],
                             "frangieh/control/gene": ["seurat_de", "liscovitch_method", "weissman_method", "seurat_de_nb"],
                             "frangieh/ifn_gamma/gene": ["seurat_de", "liscovitch_method", "weissman_method", "seurat_de_nb"],
                             "papalexi/eccite_screen/gene": ["seurat_de", "liscovitch_method", "weissman_method", "seurat_de_nb"],
                             "papalexi/eccite_screen/protein": ["seurat_de", "liscovitch_method", "weissman_method", "seurat_de_nb"],
                             "schraivogel/enhancer_screen_chr11/gene": ["seurat_de", "liscovitch_method", "weissman_method", "seurat_de_nb"],
                             "schraivogel/enhancer_screen_chr8/gene": ["seurat_de", "liscovitch_method", "weissman_method", "seurat_de_nb"]
                             // "simulated/experiment_2/gene": ["seurat_de", "liscovitch_method", "weissman_method", "seurat_de_nb"]
                            ]
*/

data_method_pairs_indiv = []

data_method_pairs_grouped = ["frangieh/co_culture/gene": ["seurat_de"],
                             "frangieh/control/gene": ["seurat_de"],
                             "frangieh/ifn_gamma/gene": ["seurat_de"],
                             "papalexi/eccite_screen/gene": ["seurat_de"],
                             "papalexi/eccite_screen/protein": ["seurat_de"],
                             "schraivogel/enhancer_screen_chr11/gene": ["seurat_de"],
                             "schraivogel/enhancer_screen_chr8/gene": ["seurat_de"],
                             "simulated/experiment_2/gene": ["seurat_de"]
                            ]

row_names = ["frangieh/co_culture/gene",
             "frangieh/control/gene",
             "frangieh/ifn_gamma/gene",
             "papalexi/eccite_screen/gene",
             "papalexi/eccite_screen/protein",
             "schraivogel/enhancer_screen_chr11/gene",
             "schraivogel/enhancer_screen_chr8/gene",
             "simulated/experiment_2/gene"]
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
                         [23,  9, 45, 26, 34, 8, 9],  // frangieh/control/gene
                         [36, 13, 13, 34, 53, 8, 13], // frangieh/ifn_gamma/gene
                         [18,  6,  6, 19, 22, 6, 6],  // papalexi/eccite_screen/gene
                         [2,   1,  1,  5,  1, 4, 1],  // papalexi/eccite_screen/protein
                         [6,  12, 12, 18,  6, 8, 12], // schraivogel/enhancer_screen_chr11/gene
                         [8,  14, 14, 21,  6, 8, 14], // schraivogel/enhancer_screen_chr8/gene
                         [8,  14, 14, 21,  6, 8, 14]  // simulated/experiment_2/gene
                         ]
// schraivogel_method, seurat_de, liscovitch_method, mimosca, weissman_method, sceptre, seurat_de_nb

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
                           ["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q"], // frangieh/co_culture/gene
                           ["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q"], // frangieh/control/gene
                           ["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q"], // frangieh/ifn_gamma/gene
                           ["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q"], // papalexi/eccite_screen/gene
                           ["short.q", "short.q", "short.q", "short.q", "short.q", "short.q", "short.q"], // papalexi/eccite_screen/protein
                           ["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q"], // schraivogel/enhancer_screen_chr11/gene
                           ["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q"], // schraivogel/enhancer_screen_chr8/gene
                           ["short.q", "short.q", "short.q", "all.q",   "short.q", "short.q", "short.q"]  // simulated/experiment_2/gene
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
