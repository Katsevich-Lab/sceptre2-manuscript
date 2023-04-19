data_method_pairs_grouped = ["frangieh/co_culture/gene": ["seurat_de", "sceptre"],
                             "frangieh/control/gene": ["seurat_de", "sceptre"],
                             "frangieh/ifn_gamma/gene": ["seurat_de"],
                             "papalexi/eccite_screen/gene": ["seurat_de"],
                             "papalexi/eccite_screen/protein": ["seurat_de", "sceptre"],
                             "schraivogel/enhancer_screen_chr11/gene": ["seurat_de", "sceptre"],
                             "schraivogel/enhancer_screen_chr8/gene": ["seurat_de", "sceptre"]
                            ]
 
data_method_pairs_indiv = []
                            
row_names = ["frangieh/co_culture/gene",
             "frangieh/control/gene",
             "frangieh/ifn_gamma/gene",
             "papalexi/eccite_screen/gene",
             "papalexi/eccite_screen/protein",
             "schraivogel/enhancer_screen_chr11/gene",
             "schraivogel/enhancer_screen_chr8/gene"]
col_names = ["seurat_de",
             "sceptre"]

// SECOND, define a matrix indicating the amount of RAM to request for each dataset-method pair
data_method_ram_matrix = [
                         [13, 8], // frangieh/co_culture/gene
                         [9,  8], // frangieh/control/gene
                         [13, 8], // frangieh/ifn_gamma/gene
                         [6,  6], // papalexi/eccite_screen/gene
                         [2,  2], // papalexi/eccite_screen/protein
                         [12, 8], // schraivogel/enhancer_screen_chr11/gene
                         [14, 8], // schraivogel/enhancer_screen_chr8/gene
                         ]
// seurat_de, sceptre

// THIRD, define a matrix indicating the queue in which to put a given dataset-method pair process
data_method_queue_matrix = [
                           ["all.q", "all.q"], // frangieh/co_culture/gene
                           ["all.q", "all.q"], // frangieh/control/gene
                           ["all.q", "all.q"], // frangieh/ifn_gamma/gene
                           ["all.q", "all.q"], // papalexi/eccite_screen/gene
                           ["all.q", "all.q"], // papalexi/eccite_screen/protein
                           ["all.q", "all.q"], // schraivogel/enhancer_screen_chr11/gene
                           ["all.q", "all.q"], // schraivogel/enhancer_screen_chr8/gene
                           ]
// seurat_de, sceptre

// FOURTH, define an ordered list of optional arguments to each of the methods (Should be strings of the form "arg1=value1;arg2=value2;arg3=value3")
optional_args = [
"", // seurat_de
"" // sceptre
]