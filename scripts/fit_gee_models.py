#!/usr/bin/env python3
"""Fit the generalised estimating equations for junction detectability.

One logistic GEE per sequencing depth, grouped by junction. The outcome is
whether a junction was detected in a given heart (DCM) patient at that depth;
the predictors are the junction and gene features held fixed at their 300M
values. Continuous predictors are z-scored on a common scale so that an odds
ratio reads as the effect of one standard deviation.

The junction universe is restricted to junctions with mean PSI between 5% and
95% at 300M that were detected in all 100 subsamples there, which is the same
set the feature analysis uses.

Takes about forty minutes and writes one row per depth and model term.

    python3 fit_gee_models.py \\
        --table comprehensive_splicing_events_mega_table_with_regions.tsv \\
        --output data/derived/gee_coefficients.tsv
"""

import argparse
import gc
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import false_discovery_control
from statsmodels.genmod.cov_struct import Independence
from statsmodels.genmod.families import Binomial

warnings.filterwarnings("ignore")

REFERENCE_DEPTH = 300
GENE_TYPE_REFERENCE = "Protein_coding"
REGION_REFERENCE = "CDS"

FORMULA = (
    "Detected ~ z_log_TPM + z_PSI_fixed + z_Compactness + z_log_gene_length "
    "+ z_GC_content + C(Gene_type) + Is_fast_evolving + Is_heart_specific "
    "+ C(Region)"
)

DTYPES = {
    "Phenotype": "str",
    "Junction": "str",
    "Samples": "uint8",
    "PSI": "float32",
    "TPM": "float32",
    "Gene_type": "str",
    "Gene_length": "float32",
    "GC_content": "float32",
    "Compactness": "float32",
    "Is_fast_evolving": "bool",
    "Is_heart_specific": "bool",
    "Region": "str",
    "Depth_numeric": "float32",
}

FEATURE_COLUMNS = ["Gene_type", "Gene_length", "GC_content", "Compactness",
                   "Is_fast_evolving", "Is_heart_specific", "Region"]


def load(path):
    print("Loading the feature table")
    raw = pd.read_csv(path, sep="\t", dtype=DTYPES, low_memory=False,
                      usecols=list(DTYPES))
    print(f"  {len(raw):,} rows, {raw['Junction'].nunique():,} junctions, "
          f"phenotypes {sorted(raw['Phenotype'].unique())}")
    return raw


def restrict(raw):
    """Junctions with PSI 5-95% and Samples == 100 at the reference depth."""
    reference = raw[raw["Depth_numeric"] == REFERENCE_DEPTH]
    mean_psi = reference.groupby("Junction")["PSI"].mean()
    within_psi = set(mean_psi[(mean_psi >= 0.05) & (mean_psi <= 0.95)].index)
    robust = set(reference[reference["Samples"] == 100]["Junction"].unique())
    keep = within_psi & robust
    print(f"  PSI 5-95%: {len(within_psi):,}   detected in 100: {len(robust):,}"
          f"   both: {len(keep):,}")
    return raw[raw["Junction"].isin(keep)].copy(), sorted(keep)


def junction_features(table):
    """One row per junction, features taken from the deepest available depth."""
    deepest_first = table.sort_values("Depth_numeric", ascending=False)
    features = (
        deepest_first.groupby("Junction", sort=False).first()[FEATURE_COLUMNS]
        .reset_index()
    )

    # PSI is held fixed at its 300M mean, falling back to the deepest depth
    # where the junction was seen.
    at_reference = (
        table[table["Depth_numeric"] == REFERENCE_DEPTH]
        .groupby("Junction")["PSI"].mean().rename("PSI_fixed").reset_index()
    )
    fallback = (
        deepest_first.groupby("Junction").first()[["PSI"]]
        .rename(columns={"PSI": "PSI_fallback"}).reset_index()
    )
    psi = at_reference.merge(fallback, on="Junction", how="right")
    psi["PSI_fixed"] = psi["PSI_fixed"].fillna(psi["PSI_fallback"]).astype("float32")

    features = features.merge(psi[["Junction", "PSI_fixed"]], on="Junction",
                              how="left")
    features["log_gene_length"] = np.log1p(
        features["Gene_length"].astype("float64")).astype("float32")
    return features


def expression(table):
    """log TPM per junction and patient, from the deepest available depth."""
    tpm = (
        table.sort_values("Depth_numeric", ascending=False)
        .groupby(["Junction", "Phenotype"], sort=False).first()[["TPM"]]
        .rename(columns={"TPM": "TPM_fixed"}).reset_index()
    )
    tpm["log_TPM"] = np.log1p(tpm["TPM_fixed"].astype("float64")).astype("float32")
    return tpm


def build_depth_frame(junctions, phenotypes, features, tpm, detected, zparams):
    """Full junction-by-patient grid for one depth, with the detection outcome."""
    grid = pd.MultiIndex.from_product(
        [junctions, phenotypes], names=["Junction", "Phenotype"]
    ).to_frame(index=False)
    grid["Detected"] = [
        int(pair in detected)
        for pair in zip(grid["Junction"], grid["Phenotype"])
    ]
    grid = grid.merge(features, on="Junction", how="left")
    grid = grid.merge(tpm[["Junction", "Phenotype", "log_TPM"]],
                      on=["Junction", "Phenotype"], how="left")
    grid["log_TPM"] = grid["log_TPM"].fillna(0.0).astype("float32")

    for column, (mean, sd) in zparams.items():
        grid[f"z_{column}"] = (
            (grid[column].astype("float64") - mean) / sd
        ).astype("float32")

    grid["Gene_type"] = pd.Categorical(
        grid["Gene_type"].fillna("Other"),
        categories=[GENE_TYPE_REFERENCE]
        + [c for c in grid["Gene_type"].unique()
           if c != GENE_TYPE_REFERENCE and not pd.isna(c)],
    )
    grid["Region"] = pd.Categorical(
        grid["Region"].fillna(REGION_REFERENCE).replace("", REGION_REFERENCE),
        categories=[REGION_REFERENCE]
        + sorted(c for c in grid["Region"].dropna().unique()
                 if c not in (REGION_REFERENCE, "")),
    )
    return grid


def coefficient_table(result, depth):
    table = result.summary2().tables[1].copy()
    table.index.name = "term"
    table = table.reset_index()
    table.columns = ["term", "coef", "std_err", "z", "p_value", "ci_low", "ci_high"]
    table["OR"] = np.exp(table["coef"])
    table["OR_low"] = np.exp(table["ci_low"])
    table["OR_high"] = np.exp(table["ci_high"])

    table["p_adj"] = np.nan
    testable = (table["term"] != "Intercept") & table["p_value"].notna()
    if testable.any():
        table.loc[testable, "p_adj"] = false_discovery_control(
            table.loc[testable, "p_value"])

    table["depth"] = depth
    table["sig"] = table["p_adj"].map(
        lambda p: "" if pd.isna(p)
        else "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
    )
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--table",
        default="comprehensive_splicing_events_mega_table_with_regions.tsv")
    parser.add_argument("--output", default="data/derived/gee_coefficients.tsv")
    args = parser.parse_args()

    raw = load(args.table)
    phenotypes = sorted(raw["Phenotype"].unique())
    depths = sorted(raw["Depth_numeric"].dropna().astype(int).unique())

    print("\nRestricting the junction universe")
    table, junctions = restrict(raw)
    del raw
    gc.collect()

    print("\nBuilding features")
    features = junction_features(table)
    tpm = expression(table)

    zparams = {
        "log_TPM": (float(tpm["log_TPM"].mean()), float(tpm["log_TPM"].std())),
        "PSI_fixed": (float(features["PSI_fixed"].mean()),
                      float(features["PSI_fixed"].std())),
        "Compactness": (float(features["Compactness"].mean()),
                        float(features["Compactness"].std())),
        "log_gene_length": (float(features["log_gene_length"].mean()),
                            float(features["log_gene_length"].std())),
        "GC_content": (float(features["GC_content"].mean()),
                       float(features["GC_content"].std())),
    }
    for column, (mean, sd) in zparams.items():
        print(f"  {column:18} mean {mean:9.4f}  sd {sd:9.4f}")

    detected_at = {}
    for depth in depths:
        rows = table[table["Depth_numeric"] == depth][["Junction", "Phenotype"]]
        detected_at[depth] = set(zip(rows["Junction"], rows["Phenotype"]))
    del table
    gc.collect()

    coefficients = {}
    for depth in [d for d in depths if d != REFERENCE_DEPTH]:
        print(f"\nDepth {depth}M")
        frame = build_depth_frame(junctions, phenotypes, features, tpm,
                                  detected_at[depth], zparams)
        print(f"  {len(frame):,} rows, detected {frame['Detected'].sum():,} "
              f"({frame['Detected'].mean():.1%})")
        model = smf.gee(FORMULA, groups="Junction", data=frame,
                        family=Binomial(), cov_struct=Independence())
        coefficients[depth] = coefficient_table(model.fit(maxiter=60), depth)
        print(coefficients[depth][["term", "OR", "p_adj", "sig"]].to_string(index=False))
        del frame
        gc.collect()

    columns = ["depth", "term", "coef", "std_err", "z", "p_value",
               "ci_low", "ci_high", "OR", "OR_low", "OR_high", "p_adj", "sig"]
    table = pd.concat(coefficients.values(), ignore_index=True)[columns]
    table.sort_values(["depth", "term"]).to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {args.output} ({len(coefficients)} depths, {len(table)} rows)")


if __name__ == "__main__":
    main()
