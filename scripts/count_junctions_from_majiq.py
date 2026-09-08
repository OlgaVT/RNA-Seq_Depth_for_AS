#!/usr/bin/env python3
"""Count the junctions involved in AS in one MAJIQ PSI table.

Prints `<file> <n_junctions> <depth> <uniquely_mapped_depth>`, the row format of
the `*_junctions` inputs in data/figure6/. Column 11 of the MAJIQ table holds
the junction coordinates of the local splicing variant, separated by semicolons.

    python3 count_junctions_from_majiq.py sample.psi.tsv 50000000 48213445
"""

import sys

majiq_path, depth, depth_unique = sys.argv[1:4]

junctions = set()
with open(majiq_path) as handle:
    for line in handle:
        fields = line.split()
        if fields[0] == "Gene":
            continue
        gene = fields[0].split(".")[0]
        for coordinates in fields[10].split(";"):
            junctions.add(gene + "_" + coordinates)

print(majiq_path, len(junctions), depth, depth_unique)
