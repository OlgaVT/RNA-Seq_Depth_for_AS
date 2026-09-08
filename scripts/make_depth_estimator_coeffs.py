#!/usr/bin/env python3
"""
make_coeffs.py

Precompute coefficients for sequencing depth estimation from:
  Input_Files/J26675_summary.tsv       (dilated cardiomyopathy)
  Input_Files/hypothalamus_summary.tsv (hypothalamus)
  Input_Files/adipose_1_summary.tsv    (adipose pre-treatment)
  Input_Files/adipose_2_summary.tsv    (adipose post-treatment)

and write a compact JSON for use in a browser:

data/seq_depth_coeffs.json

JSON structure (roughly):

{
  "meta": { ... },
  "data": {
    "dilated cardiomyopathy": {
      "Gene": {
        "0": { depths, counts, coef1, coef2, x0 },
        "1": { ... },
        ...
      },
      "Junction": {
        ...
      }
    },
    "hypothalamus": { ... },
    "adipose pre-treatment": { ... },
    "adipose post-treatment": { ... }
  }
}
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------- CONFIG ----------------

INPUT_DIR = Path("Input_Files")

# Mapping from human-readable dataset name (used in the web UI)
# to the *file stem* (without "_summary.tsv").
DATASETS = {
    "dilated cardiomyopathy": "J26675",
    "hypothalamus": "hypothalamus",
    "adipose pre-treatment": "adipose_1",
    "adipose post-treatment": "adipose_2",
}

# Where to write JSON (relative to current working directory)
OUTPUT_JSON = Path("data/seq_depth_coeffs.json")

# Confidences to precompute
CONF_RANGE = range(0, 101)  # 0..100 inclusive

# Types we support; must match column names in your TSVs
TYPES = ("Gene", "Junction")


# ------------- HELPERS -------------------

def safe_polyfit(x, y, deg):
    """
    Polynomial fit with basic safety.
    Returns list of floats (len deg+1) or [None,...] if fit fails.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < (deg + 1):
        return [None] * (deg + 1)
    try:
        coefs = np.polyfit(x[mask], y[mask], deg)
        return [float(c) for c in coefs]
    except Exception:
        return [None] * (deg + 1)


def read_summary_file(path: Path) -> pd.DataFrame:
    """
    Read a summary TSV: whitespace-delimited with columns like:
    Phenotype  Gene  Junction  Depth  Samples  PSI  TPM
    """
    # Your files are whitespace-separated
    df = pd.read_csv(path, sep=r"\s+", engine="python")

    # Normalize column names a bit
    df.columns = df.columns.str.strip()

    required = {"Depth", "Samples"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {missing}. Found columns: {list(df.columns)}")

    # Samples numeric
    df["Samples"] = pd.to_numeric(df["Samples"], errors="coerce")
    df = df.dropna(subset=["Samples"])

    return df


def aggregate_counts(df: pd.DataFrame, confidence: int, type_col: str):
    """
    For a given dataframe, confidence and type ('Gene' or 'Junction'),
    implements the logic:

      filtered_df = df[df["Samples"] >= confidence]
      filtered_df = filtered_df.drop_duplicates(subset=["Depth", type_col])
      grouped     = filtered_df.groupby("Depth")[type_col].count()

    Then extract:

      depths: numeric part of Depth strings, sorted ascending (e.g. [50,100,150,...])
      counts: counts at each depth, ordered to match depths
      x_norm: counts / max_count

    And compute linear & quadratic polyfits y(depth_in_M) = f(x_norm),
    where x_norm ~ fraction of events detected.
    """
    if type_col not in df.columns:
        return None

    filtered = df[df["Samples"] >= confidence].copy()
    if filtered.empty:
        return None

    filtered = filtered.drop_duplicates(subset=["Depth", type_col])

    grouped = filtered.groupby("Depth")[type_col].count()
    if grouped.empty:
        return None

    tmp = grouped.reset_index()  # columns: Depth, type_col
    tmp["depth_M"] = tmp["Depth"].astype(str).str.extract(r"(\d+)").astype(float)
    tmp = tmp.dropna(subset=["depth_M"]).sort_values("depth_M")
    if tmp.empty:
        return None

    depths = tmp["depth_M"].tolist()
    counts = tmp[type_col].astype(int).tolist()

    max_count = counts[-1] if counts else 0
    if max_count <= 0:
        return None

    x_norm = [c / max_count for c in counts]
    y = depths

    coef1 = safe_polyfit(x_norm, y, 1)  # linear
    coef2 = safe_polyfit(x_norm, y, 2)  # quadratic

    if all(c is None for c in (coef1 + coef2)):
        return None

    scenario = {
        "depths": depths,   # [50,100,150,...] in M reads
        "counts": counts,   # counts at each depth
        "coef1": coef1,     # [a1, b1] for y = a1 * p + b1
        "coef2": coef2,     # [a2, b2, c2] for y = a2 * p^2 + b2 * p + c2
        "x0": x_norm[0],    # normalized fraction at the first depth
    }
    return scenario


def build_dataset_entry(human_name: str, filestem: str):
    """
    For a given dataset label (human_name) and file stem,
    read Input_Files/{filestem}_summary.tsv and compute per-type coeffs.
    Returns a dict: { "Gene": {conf->scenario}, "Junction": {conf->scenario}, ... }
    """
    path = INPUT_DIR / f"{filestem}_summary.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Expected file not found: {path}")

    df = read_summary_file(path)
    entry = {}

    for type_col in TYPES:
        if type_col not in df.columns:
            continue

        conf_map = {}
        for conf in CONF_RANGE:
            scenario = aggregate_counts(df, conf, type_col)
            if scenario is not None:
                conf_map[str(conf)] = scenario

        if conf_map:
            entry[type_col] = conf_map

    return entry


def main():
    out_data = {}

    for human_name, stem in DATASETS.items():
        print(f"Processing dataset: {human_name} (file stem: {stem})")
        entry = build_dataset_entry(human_name, stem)
        if entry:
            out_data[human_name] = entry
        else:
            print(f"  Warning: no usable data for {human_name}")

    if not out_data:
        raise SystemExit("No usable datasets. Nothing to write.")

    out = {
        "meta": {
            "datasets": list(out_data.keys()),
            "types": list(TYPES),
            "confidence_min": min(CONF_RANGE),
            "confidence_max": max(CONF_RANGE),
            "description": (
                "For each dataset (human-readable name) and type (Gene/Junction), "
                "stores, for each confidence, depths (M reads), counts, and polyfit "
                "coefficients for read depth estimation."
            ),
        },
        "data": out_data,
    }

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_JSON, "w") as f:
        json.dump(out, f, indent=2)

    print(f"Wrote {OUTPUT_JSON.resolve()}")


if __name__ == "__main__":
    main()
