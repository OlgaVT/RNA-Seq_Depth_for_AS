# Inputs

Nothing here except `derived/` is written by this repository.

## expression_category_totals.tsv

Genes and junctions per TPM category and cohort, with the standard deviation
over samples where there are several. It comes from the featureCounts TPM
matrices and the regtools junction lists, upstream of the splicing analysis and
too large to keep here.

It counts every expressed gene in a category, not only those with alternative
splicing, which is why it cannot come out of the MAJIQ summaries. The notebooks
use it as the denominator for the increment tables.

## annotation/

Per-gene properties the junction feature table is built from: type, length,
compactness and GC content from Ensembl and GENCODE, the genes overlapping
human accelerated regions from `scripts/har_genes.sh`, and the heart-elevated
gene list.

## derived/

| File | Built by |
|---|---|
| `feature_counts.tsv`, `feature_associations*.tsv` | `scripts/summarise_junction_features.py` |
| `gee_coefficients.tsv` | `scripts/fit_gee_models.py` |
| `extrapolation_counts/` | `scripts/count_genes_from_majiq.py`, `count_junctions_from_majiq.py`, one line per MAJIQ output file |
| `cost_counts/` | detections per subsample count and expression level, from the MAJIQ output |
| `nease/` | NEASE, run per sample and depth on the DICAST-unified events |
