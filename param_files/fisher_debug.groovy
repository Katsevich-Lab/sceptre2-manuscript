data_method_pairs = ["schraivogel/ground_truth_perturbseq/gene": ["fisher_exact"],
                     "schraivogel/ground_truth_tapseq/gene": ["fisher_exact"]
                     ]

row_names = ["schraivogel/ground_truth_perturbseq/gene",
             "schraivogel/ground_truth_tapseq/gene"]
col_names = ["fisher_exact_thresholded"]

data_method_ram_matrix = [
[5], // schraivogel/ground_truth_perturbseq/gene
[5], // schraivogel/ground_truth_tapseq/gene
]

data_method_queue_matrix = [
["short.q"], // schraivogel/ground_truth_perturbseq/gene
["short.q"], // schraivogel/ground_truth_tapseq/gene
]

optional_args = [
"" // fisher_exact
]
