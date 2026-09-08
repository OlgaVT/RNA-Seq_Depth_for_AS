#!/usr/bin/env python3
"""Normalise a featureCounts output to TPM.

Writes `<title>_TPM` next to the working directory.

    python3 tpm.py featurecounts_output sample_title
"""

import sys

import pandas as pd

counts_path, title = sys.argv[1:3]

table = pd.read_csv(counts_path, sep="\t", skiprows=1)
table.columns = ["Geneid", "Chr", "Start", "End", "Strand", "Length", title]

reads_per_kilobase = table[title] * 1000 / table["Length"]
table["TPM"] = reads_per_kilobase / (reads_per_kilobase.sum() / 1_000_000)

table.to_csv(f"{title}_TPM", sep="\t", index=False)
