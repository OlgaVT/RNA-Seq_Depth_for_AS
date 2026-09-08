#!/usr/bin/env python3
"""Count the genes with AS in one MAJIQ PSI table.

Prints `<file> <n_genes> <depth> <uniquely_mapped_depth>`, the row format of the
`*_genes` inputs in data/figure6/.

    python3 count_genes_from_majiq.py sample.psi.tsv 50000000 48213445
"""

import sys

majiq_path, depth, depth_unique = sys.argv[1:4]

genes = set()
with open(majiq_path) as handle:
    for line in handle:
        first = line.split()[0]
        if first != "Gene":
            genes.add(first.split(".")[0])

print(majiq_path, len(genes), depth, depth_unique)
