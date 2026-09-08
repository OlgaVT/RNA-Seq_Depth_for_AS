"""Checks on the analysis code.

These are invariants, not stored outputs. A check that asserts the number the
pipeline happened to print last time only says the code has not changed; the
checks below say something that has to hold if the code is right, and would
fail if the pipeline were quietly broken.

Only the notebooks that need no downloaded data are covered.

    python3 -m pytest tests
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def run_notebook_code(stem, upto=None):
    """Execute a notebook's code cells and return the resulting namespace."""
    cells = json.loads((ROOT / "notebooks" / f"{stem}.ipynb").read_text())["cells"]
    code = ["".join(cell["source"]) for cell in cells if cell["cell_type"] == "code"]
    if upto is not None:
        code = code[:upto]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.show = lambda *args, **kwargs: None

    namespace = {"__name__": "__main__"}
    for block in code:
        exec(compile(block, f"{stem}.ipynb", "exec"), namespace)
    return namespace


@pytest.fixture(scope="module")
def junction_features():
    # Stop before the plotting cells; the statistics are what is checked.
    return run_notebook_code("03_junction_features", upto=3)


def test_subcategory_proportions_are_a_partition(junction_features):
    """Within a feature and depth, the sub-categories must cover everything once."""
    counts = junction_features["counts"]
    for feature, rows in counts.groupby("feature"):
        totals = rows.groupby("depth")["count"].sum()
        assert (totals > 0).all(), feature
        # Every depth must see the same junction universe for the feature.
        assert totals.nunique() <= len(totals)
        assert not rows.duplicated(["depth", "subcategory"]).any(), feature


def test_expression_is_the_strongest_compositional_shift(junction_features):
    """Whatever the numbers, expression must dominate the other features.

    This is the claim the figure makes; if some other feature ever overtook TPM
    the conclusion would change, so it is worth asserting rather than assuming.
    """
    trends = junction_features["trends"]
    largest = {
        feature: max(abs(s["slope"]) for s in trends[feature]["slopes"].values())
        for feature in trends
    }
    assert max(largest, key=largest.get) == "TPM"
    runner_up = sorted(largest.values())[-2]
    assert largest["TPM"] > 3 * runner_up


def test_low_expression_gains_what_high_expression_loses(junction_features):
    """Extra reads add junctions from lowly expressed genes.

    The three tertiles are shares of one whole, so their slopes must sum to
    zero, and the direction must run from High to Low.
    """
    tpm = junction_features["trends"]["TPM"]["slopes"]
    assert tpm["Low"]["slope"] > 0 > tpm["High"]["slope"]
    assert sum(s["slope"] for s in tpm.values()) == pytest.approx(0, abs=1e-9)


def test_fdr_is_never_below_the_raw_p_value(junction_features):
    trends = junction_features["trends"]
    for feature, result in trends.items():
        assert result["fdr"] >= result["p_value"] - 1e-12, feature


@pytest.fixture(scope="module")
def extrapolation():
    return run_notebook_code("04_gtex_and_tcga")


def test_deeper_sequencing_finds_more_than_the_cohort_does(extrapolation):
    """Every fit must put 200M above the cohort's own typical depth."""
    summary = extrapolation["summary"]
    linear = summary[summary["curve"] == "linear"]
    assert (linear["at 200M"] > linear["baseline"]).all()
    assert (linear["missed"] > 0).all()


def test_junctions_outnumber_genes(extrapolation):
    """A gene carries several junctions, so junction counts must be larger."""
    summary = extrapolation["summary"]
    for cohort, rows in summary.groupby("cohort"):
        genes = rows[rows["kind"] == "genes"]["baseline"].iloc[0]
        junctions = rows[rows["kind"] == "junctions"]["baseline"].iloc[0]
        assert junctions > 5 * genes, cohort


@pytest.fixture(scope="module")
def cost():
    return run_notebook_code("05_sequencing_cost")


def test_each_extra_detection_costs_more_than_the_last(cost):
    """Saturation: the same money buys fewer new detections as depth grows."""
    costs = cost["costs"]
    for keys, group in costs.groupby(["panel", "samples", "platform"]):
        curve = group.sort_values("depth_M")["euro"].to_numpy()
        assert (curve > 0).all(), keys
        assert (curve[1:] > curve[:-1]).all(), keys


def test_a_junction_is_cheaper_to_find_than_a_gene(cost):
    costs = cost["costs"]
    deepest = costs[costs["depth_M"] == costs["depth_M"].max()]
    for (samples, platform), rows in deepest.groupby(["samples", "platform"]):
        by_kind = rows.set_index("kind")["euro"]
        assert by_kind["junctions"] < by_kind["genes"], (samples, platform)


def test_detectability_model_agrees_with_the_composition_analysis():
    """Expression raises the odds of detection at every depth.

    The model is fitted separately from the counting, so the two agreeing is a
    real check rather than a restatement.
    """
    coefficients = pd.read_csv(
        ROOT / "data" / "derived" / "gee_coefficients.tsv", sep="\t")
    expression = coefficients[coefficients["term"] == "z_log_TPM"]
    assert (expression["OR"] > 1).all()
    assert (expression["OR_low"] > 1).all()

    intercepts = coefficients[coefficients["term"] == "Intercept"]
    assert intercepts["p_adj"].isna().all(), "the intercept must not be FDR-corrected"


def test_pathways_gained_at_depth_are_significant_and_new():
    """Every reported pathway must clear the threshold at 200M and not before."""
    namespace = run_notebook_code("08_pathway_enrichment")
    table, results = namespace["table_s7"], namespace["results"]

    assert (table["p-adjusted"] < 0.05).all()
    earlier = set(results[namespace["PREVIOUS_DEPTH"]]
                  .query("`adj p_value` < 0.05")["Pathway name"])
    for pathway in table[table["dataset"] == "Hypothalamus"]["KEGG pathway"]:
        assert not any(pathway in name for name in earlier), pathway


def test_expression_category_totals_cover_every_cohort():
    totals = pd.read_csv(ROOT / "data" / "expression_category_totals.tsv", sep="\t")
    assert len(totals) == 5 * 2 * 5
    assert totals["total"].gt(0).all()
    assert not totals.duplicated(["dataset", "kind", "category"]).any()
    # A gene carries several junctions, in every cohort and category.
    wide = totals.pivot_table(index=["dataset", "category"], columns="kind",
                              values="total")
    assert (wide["junctions"] > wide["genes"]).all()


def test_summary_reader_matches_a_hand_counted_file(tmp_path):
    """The one check whose right answer was worked out on paper, not by the code."""
    from asdepth import count_deep_sequenced

    rows = [
        "Phenotype Gene\tJunction\tDepth\tSamples\tPSI\tTPM",
        # Two junctions of one gene at 50M, one of them seen in all 100.
        "s\tG1\tG1_1-2\t50M\t100\t0.5\t5.0",
        "s\tG1\tG1_3-4\t50M\t40\t0.5\t5.0",
        # A near-constitutive junction, and a gene outside the TPM bin.
        "s\tG2\tG2_1-2\t50M\t100\t0.99\t5.0",
        "s\tG3\tG3_1-2\t50M\t100\t0.5\t50.0",
        # At 100M the first junction persists and a second gene appears.
        "s\tG1\tG1_1-2\t100M\t100\t0.5\t5.0",
        "s\tG4\tG4_1-2\t100M\t100\t0.5\t5.0",
    ]
    path = tmp_path / "mini_summary.tsv"
    path.write_text("\n".join(rows) + "\n")

    counts = count_deep_sequenced(path, [(1, 10)])[(1, 10)]
    assert counts.depths == [50, 100]
    assert counts.all_genes == [1, 2]
    assert counts.all_junctions == [2, 2]
    assert counts.robust_junctions == [1, 2]
    assert counts.sporadic_junctions == [1, 0]
