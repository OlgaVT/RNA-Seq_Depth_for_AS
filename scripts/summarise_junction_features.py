#!/usr/bin/env python3
"""Reduce the junction feature table to the summaries the notebooks plot.

Two summaries come out of it, both small enough to keep in the repository:

    data/derived/feature_counts.tsv
        junctions per sequencing depth and feature sub-category. Continuous
        features are split into tertiles, categorical ones kept as they are.
    data/derived/feature_associations*.tsv
        a mixed association matrix over the nine features, at the deepest
        sequencing depth: signed Pearson r for numeric pairs, the correlation
        ratio for continuous against nominal, Cramer's V for nominal pairs.

Only junctions detected in all 100 subsamples are counted, so that subsampling
noise does not drive the composition.

    python3 scripts/summarise_junction_features.py \\
        --table comprehensive_splicing_events_mega_table_with_regions.tsv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import chi2_contingency

CONTINUOUS = ["PSI", "GC_content", "Compactness", "TPM"]
CATEGORICAL = ["Gene_type", "Region", "Is_heart_specific", "Is_fast_evolving"]
FEATURES = CONTINUOUS + CATEGORICAL

COLUMNS = ["Depth", "Samples", "PSI", "TPM", "Gene_type", "Region",
           "GC_content", "Compactness", "Is_fast_evolving", "Is_heart_specific"]
DTYPES = {"Depth": "category", "Samples": "uint8", "PSI": "float32",
          "TPM": "float32", "GC_content": "float32", "Compactness": "float32"}

# The association matrix carries two derived columns beyond the eight features:
# log TPM, and how far PSI sits from an even split.
PANEL = ["log_TPM", "PSI", "PSI_ext", "Gene_type", "Compactness", "Region",
         "Is_fast_evolving", "Is_heart_specific", "GC_content"]
KINDS = {"log_TPM": "continuous", "PSI": "continuous", "PSI_ext": "continuous",
         "Compactness": "continuous", "GC_content": "continuous",
         "Gene_type": "nominal", "Region": "nominal",
         "Is_fast_evolving": "boolean", "Is_heart_specific": "boolean"}
ABBREVIATIONS = {"Protein_coding": "PC", "lncRNA": "lnc", "Pseudogene": "Ps",
                 "Other": "Oth", "CDS": "CDS", "5UTR": "5'U", "3UTR": "3'U"}


def load(path):
    table = pd.read_csv(path, sep="\t", usecols=COLUMNS, dtype=DTYPES)
    robust = table[table["Samples"] == 100].copy()
    robust["depth"] = robust["Depth"].astype(str).str[:-1].astype(int)
    robust = robust.drop(columns=["Depth", "Samples"])
    robust = robust.dropna(how="all", subset=FEATURES)
    print(f"{len(robust):,} junction observations, "
          f"depths {sorted(robust['depth'].unique())}")
    return robust


def subcategory_counts(table, feature, continuous):
    """Junctions per sequencing depth and sub-category of one feature."""
    data = table[table[feature].notna()]
    if continuous:
        subcategory = pd.qcut(
            data[feature], q=3, labels=["Low", "Medium", "High"], duplicates="drop"
        )
    else:
        subcategory = data[feature]

    counts = (
        data.groupby(["depth", subcategory.rename("subcategory")], observed=True)
        .size().rename("count").reset_index()
    )
    counts.insert(0, "feature", feature)
    return counts


def correlation_ratio(categories, values):
    mask = ~(pd.isnull(categories) | np.isnan(values))
    categories, values = categories[mask], values[mask]
    grand_mean = values.mean()
    between = sum(
        ((values[categories == c]).mean() - grand_mean) ** 2 * (categories == c).sum()
        for c in np.unique(categories)
    )
    total = ((values - grand_mean) ** 2).sum()
    return np.sqrt(between / total) if total > 0 else 0.0


def cramers_v(a, b):
    mask = ~(pd.isnull(a) | pd.isnull(b))
    table = pd.crosstab(a[mask], b[mask])
    chi2 = chi2_contingency(table)[0]
    minimum = min(table.shape) - 1
    n = table.values.sum()
    return np.sqrt(chi2 / (n * minimum)) if n > 0 and minimum > 0 else 0.0


def association(frame, first, second):
    """Association value, and whether it carries a sign."""
    kinds = {KINDS[first], KINDS[second]}
    a, b = frame[first].values, frame[second].values
    if kinds <= {"continuous", "boolean"}:
        mask = ~(np.isnan(a.astype(float)) | np.isnan(b.astype(float)))
        return stats.pearsonr(a[mask].astype(float), b[mask].astype(float))[0], True
    if kinds == {"nominal"}:
        return cramers_v(a, b), False
    if "nominal" in kinds and "continuous" in kinds:
        values, categories = (a, b) if KINDS[first] == "continuous" else (b, a)
        return correlation_ratio(categories, values.astype(float)), False
    return cramers_v(a.astype(str), b.astype(str)), False


def dominant_category(frame, continuous, nominal):
    """The nominal class with the highest mean of the continuous feature."""
    values = frame[continuous].astype(float).values
    categories = frame[nominal].astype(object).values
    mask = ~(np.isnan(values) | pd.isnull(categories))
    values, categories = values[mask], categories[mask]
    means = {c: values[categories == c].mean() for c in np.unique(categories)}
    top = max(means, key=means.get)
    return ABBREVIATIONS.get(str(top), str(top)[:4])


def associations(table):
    deepest = table[table["depth"] == table["depth"].max()].copy()
    deepest["log_TPM"] = np.log1p(deepest["TPM"].astype(float))
    deepest["PSI_ext"] = (0.5 - deepest["PSI"].astype(float)).abs()

    size = len(PANEL)
    matrix = np.eye(size)
    signed = np.ones((size, size), dtype=bool)
    labels = np.full((size, size), "", dtype=object)

    for i, first in enumerate(PANEL):
        for j in range(i + 1, size):
            second = PANEL[j]
            value, is_signed = association(deepest, first, second)
            matrix[i, j] = matrix[j, i] = value
            signed[i, j] = signed[j, i] = is_signed
            if (KINDS[first] == "nominal") != (KINDS[second] == "nominal"):
                continuous = first if KINDS[first] == "continuous" else second
                nominal = second if KINDS[first] == "continuous" else first
                if KINDS[continuous] == "continuous":
                    label = dominant_category(deepest, continuous, nominal)
                    labels[i, j] = labels[j, i] = label

    frame = lambda values: pd.DataFrame(values, index=PANEL, columns=PANEL)
    return frame(matrix), frame(signed), frame(labels), int(deepest["depth"].iloc[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--table",
        default="comprehensive_splicing_events_mega_table_with_regions.tsv")
    parser.add_argument("--out-dir", default="data/derived")
    args = parser.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    table = load(args.table)

    counts = pd.concat(
        [subcategory_counts(table, f, f in CONTINUOUS) for f in FEATURES],
        ignore_index=True,
    )
    counts.to_csv(out / "feature_counts.tsv", sep="\t", index=False)
    print(f"wrote {out / 'feature_counts.tsv'}: {len(counts)} rows")

    matrix, signed, labels, depth = associations(table)
    matrix.to_csv(out / "feature_associations.tsv", sep="\t")
    signed.to_csv(out / "feature_associations_signed.tsv", sep="\t")
    labels.to_csv(out / "feature_associations_labels.tsv", sep="\t")
    print(f"wrote the association matrix at {depth}M reads")


if __name__ == "__main__":
    main()
