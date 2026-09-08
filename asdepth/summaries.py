"""Readers for the MAJIQ summary tables.

Two layouts are used in the study.

Deep-sequenced cohorts (adipose, heart, hypothalamus), one row per junction per
downsampling depth, tab separated with a header line::

    Phenotype Gene<TAB>Junction<TAB>Depth<TAB>Samples<TAB>PSI<TAB>TPM

`Samples` counts in how many of the 100 subsamples at that depth the junction
was detected.

SARS-CoV-2 cohort, one row per junction per sample, no header::

    Phenotype<TAB>Sample<TAB>Gene<TAB>Junction<TAB>PSI<TAB>Depth<TAB>TPM

The files run to 2 GB, so they are streamed line by line and all expression
categories are counted in a single pass.
"""

from collections import defaultdict
from dataclasses import dataclass, field

PSI_RANGE = (0.05, 0.95)


@dataclass
class DepthCounts:
    """Detections per sequencing depth for one expression category.

    `all_*` counts genes or junctions seen in at least one of the 100
    subsamples, `robust_*` only those seen in all 100, and `sporadic_*` is the
    difference, which the yellow segment of the stacked bars shows.
    """

    depths: list = field(default_factory=list)
    all_genes: list = field(default_factory=list)
    robust_genes: list = field(default_factory=list)
    sporadic_genes: list = field(default_factory=list)
    all_junctions: list = field(default_factory=list)
    robust_junctions: list = field(default_factory=list)
    sporadic_junctions: list = field(default_factory=list)


def category_of(value, categories):
    """First half-open [low, high) interval in `categories` containing `value`."""
    for low, high in categories:
        if low <= value < high:
            return low, high
    return None


def count_deep_sequenced(path, categories, psi_range=PSI_RANGE, genes=None):
    """Count genes and junctions with AS per depth, for every expression category.

    `genes` optionally restricts the count to a set of Ensembl gene IDs.
    Returns a dict mapping each category to a `DepthCounts`.
    """
    psi_low, psi_high = psi_range
    all_genes = {c: defaultdict(set) for c in categories}
    all_junctions = {c: defaultdict(set) for c in categories}
    robust_genes = {c: defaultdict(set) for c in categories}
    robust_junctions = {c: defaultdict(set) for c in categories}

    with open(path) as handle:
        for line in handle:
            phenotype, gene, junction, depth, samples, psi, tpm = line.split()
            if phenotype == "Phenotype":
                continue
            if genes is not None and gene not in genes:
                continue
            if not psi_low <= float(psi) < psi_high:
                continue
            category = category_of(float(tpm), categories)
            if category is None:
                continue

            depth = int(depth[:-1])
            all_genes[category][depth].add(gene)
            all_junctions[category][depth].add(junction)
            if int(samples) == 100:
                robust_genes[category][depth].add(gene)
                robust_junctions[category][depth].add(junction)

    results = {}
    for category in categories:
        counts = DepthCounts()
        for depth in sorted(all_genes[category]):
            genes_all = all_genes[category][depth]
            genes_robust = robust_genes[category][depth]
            junctions_all = all_junctions[category][depth]
            junctions_robust = robust_junctions[category][depth]
            counts.depths.append(depth)
            counts.all_genes.append(len(genes_all))
            counts.robust_genes.append(len(genes_robust))
            counts.sporadic_genes.append(len(genes_all - genes_robust))
            counts.all_junctions.append(len(junctions_all))
            counts.robust_junctions.append(len(junctions_robust))
            counts.sporadic_junctions.append(len(junctions_all - junctions_robust))
        results[category] = counts
    return results


def count_per_sample(path, categories, psi_range=PSI_RANGE):
    """Per-sample gene and junction counts for the SARS-CoV-2 cohort.

    Returns `(counts, depth_of_sample)`. `counts[category]['genes'][sample]` is
    the number of genes with AS in that sample; `depth_of_sample` maps a sample
    name to its depth bin ('60M' ... '150M').
    """
    psi_low, psi_high = psi_range
    genes = {c: defaultdict(set) for c in categories}
    junctions = {c: defaultdict(set) for c in categories}
    depth_of_sample = {}

    with open(path) as handle:
        for line in handle:
            phenotype, sample, gene, junction, psi, depth, tpm = line.split()
            depth_of_sample[sample] = depth
            if not psi_low <= float(psi) < psi_high:
                continue
            category = category_of(float(tpm), categories)
            if category is None:
                continue
            genes[category][sample].add(gene)
            junctions[category][sample].add(gene + "_" + junction)

    counts = {
        category: {
            "genes": {s: len(v) for s, v in genes[category].items()},
            "junctions": {s: len(v) for s, v in junctions[category].items()},
        }
        for category in categories
    }
    return counts, depth_of_sample


def read_depth_counts(path):
    """Read a `<name>_genes` / `<name>_junctions` count file.

    Space separated, one row per MAJIQ output file:
    `<file> <count> <depth> <uniquely_mapped_depth>`.
    """
    import pandas as pd

    table = pd.read_csv(path, sep=" ", header=None, index_col=0)
    table.columns = ["count", "depth", "depth_unique"]
    return table
