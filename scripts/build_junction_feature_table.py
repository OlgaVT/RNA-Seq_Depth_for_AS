#!/usr/bin/env python3
"""Build the junction feature table behind the feature and detectability analyses.

One row per junction, heart (DCM) sample and downsampling depth, annotated with
the properties of the junction and of its gene.

Inputs
    <input-dir>/J2667{5,6,7,8}_summary.tsv   MAJIQ summaries
    gencode.v44.annotation.gtf               CDS and UTR coordinates
    data/annotation/                         gene types, gene length,
                                             compactness, GC content, the genes
                                             overlapping HARs and the
                                             heart-elevated genes

The output is roughly 1 GB, which is why it is not in the repository.

    python3 build_junction_feature_table.py --input-dir Input_Files \\
        --annotation gencode.v44.annotation.gtf \\
        --output comprehensive_splicing_events_mega_table_with_regions.tsv
"""

import argparse
import re
from collections import defaultdict

import numpy as np
import pandas as pd

HEART_SAMPLES = ["J26675", "J26676", "J26677", "J26678"]

GENE_ID = re.compile(r'gene_id "([^"]+)"')
TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')

# A junction is assigned to the region it overlaps most; if it overlaps none,
# a boundary within this many bases of a region still counts.
BOUNDARY_TOLERANCE = 100


def parse_regions(gtf_path):
    """CDS, 5'UTR and 3'UTR spans per gene.

    GENCODE labels both UTRs simply 'UTR', so they are separated by their
    position relative to the coding span, taking the strand into account.
    """
    per_gene = defaultdict(lambda: {"CDS": [], "UTR": [], "strand": None})

    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("CDS", "UTR"):
                continue

            gene = GENE_ID.search(fields[8])
            transcript = TRANSCRIPT_ID.search(fields[8])
            if not gene or not transcript:
                continue

            entry = per_gene[gene.group(1).split(".")[0]]
            entry[fields[2]].append((int(fields[3]), int(fields[4])))
            entry["strand"] = fields[6]

    regions = {}
    for gene, entry in per_gene.items():
        if not entry["CDS"]:
            continue
        cds_start = min(start for start, _ in entry["CDS"])
        cds_end = max(end for _, end in entry["CDS"])
        spans = {"CDS": (cds_start, cds_end)}

        upstream, downstream = [], []
        for start, end in entry["UTR"]:
            if end <= cds_start:
                (upstream if entry["strand"] == "+" else downstream).append((start, end))
            elif start >= cds_end:
                (downstream if entry["strand"] == "+" else upstream).append((start, end))

        for name, utrs in (("5UTR", upstream), ("3UTR", downstream)):
            if utrs:
                spans[name] = (min(s for s, _ in utrs), max(e for _, e in utrs))
        regions[gene] = spans
    return regions


def junction_region(junction, gene, regions):
    """Region a junction falls in, or None if the gene has no annotation."""
    spans = regions.get(gene)
    if spans is None:
        return None
    try:
        start, end = (int(value) for value in junction.split("_")[1].split("-"))
    except (IndexError, ValueError):
        return None

    best, best_overlap, fallback = None, 0, None
    for name in ("5UTR", "CDS", "3UTR"):
        span = spans.get(name)
        if span is None:
            continue
        overlap = min(end, span[1]) - max(start, span[0])
        if overlap > best_overlap:
            best, best_overlap = name, overlap
        if fallback is None and (
            span[0] - BOUNDARY_TOLERANCE <= start <= span[1] + BOUNDARY_TOLERANCE
            or span[0] - BOUNDARY_TOLERANCE <= end <= span[1] + BOUNDARY_TOLERANCE
        ):
            fallback = name
    return best or fallback


def gene_type_group(value):
    """Collapse the Ensembl gene biotypes into the five categories used."""
    if pd.isna(value):
        return "Other"
    value = str(value).lower()
    if "protein_coding" in value or "protein coding" in value:
        return "Protein_coding"
    if any(token in value for token in ("lncrna", "lincrna", "long_ncrna")):
        return "lncRNA"
    if "pseudo" in value:
        return "Pseudogene"
    if any(token in value for token in
           ("mirna", "scrna", "snrna", "snorna", "rrna", "trna", "ncrna")):
        return "Small_RNA"
    return "Other"


def read_gene_set(path):
    with open(path) as handle:
        return {line.strip().split(".")[0] for line in handle if line.strip()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="Input_Files")
    parser.add_argument("--annotation", default="gencode.v44.annotation.gtf")
    parser.add_argument("--gene-types", default="data/annotation/gene_types.tsv")
    parser.add_argument("--gene-stats", default="data/annotation/gene_stats.tsv")
    parser.add_argument("--har-genes", default="data/annotation/har_genes.txt")
    parser.add_argument("--heart-genes",
                        default="data/annotation/heart_specific_genes.txt")
    parser.add_argument(
        "--output",
        default="comprehensive_splicing_events_mega_table_with_regions.tsv",
    )
    args = parser.parse_args()

    print("Reading summaries")
    frames = []
    for sample in HEART_SAMPLES:
        frame = pd.read_csv(f"{args.input_dir}/{sample}_summary.tsv", sep="\t")
        frame = frame.reset_index()
        frame.columns = ["Phenotype", "Gene", "Junction", "Depth",
                         "Samples", "PSI", "TPM"]
        frames.append(frame)
        print(f"  {sample}: {len(frame):,} rows")
    table = pd.concat(frames, ignore_index=True)
    table["Gene_clean"] = table["Gene"].str.split(".").str[0]

    print("Parsing the annotation")
    regions = parse_regions(args.annotation)
    print(f"  {len(regions):,} genes with a coding span")

    # Annotate one junction at a time rather than one row at a time; the same
    # junction appears at every depth and in every sample.
    unique = table[["Junction", "Gene_clean"]].drop_duplicates()
    region_of = {
        junction: junction_region(junction, gene, regions)
        for junction, gene in zip(unique["Junction"], unique["Gene_clean"])
    }
    table["Region"] = table["Junction"].map(region_of)

    print("Adding gene features")
    gene_types = pd.read_csv(args.gene_types, sep="\t")
    biotype = dict(zip(gene_types["Gene stable ID"], gene_types["Gene type"]))
    table["Gene_type_raw"] = table["Gene_clean"].map(biotype)
    table["Gene_type"] = table["Gene_type_raw"].map(gene_type_group)

    stats = pd.read_csv(args.gene_stats, sep="\t").set_index("gene_id")
    for column, source in (
        ("GC_content", "Gene % GC content"),
        ("Compactness", "compactness"),
        ("Gene_length", "gene_length"),
    ):
        table[column] = table["Gene_clean"].map(stats[source])

    table["Is_fast_evolving"] = table["Gene_clean"].isin(read_gene_set(args.har_genes))
    table["Is_heart_specific"] = table["Gene_clean"].isin(
        read_gene_set(args.heart_genes))
    table["Depth_numeric"] = table["Depth"].str.rstrip("M").astype(int)

    columns = [
        "Phenotype", "Gene", "Gene_clean", "Junction", "Depth", "Depth_numeric",
        "Samples", "PSI", "TPM", "Gene_type", "Gene_type_raw", "Region",
        "GC_content", "Compactness", "Gene_length",
        "Is_fast_evolving", "Is_heart_specific",
    ]
    table[columns].to_csv(args.output, sep="\t", index=False)

    print(f"\nWrote {args.output}: {len(table):,} rows, "
          f"{table['Junction'].nunique():,} junctions, "
          f"{table['Gene_clean'].nunique():,} genes")
    print(table["Region"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
