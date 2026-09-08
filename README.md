# RNA-Seq_Depth_for_AS

How deep does an RNA-Seq library have to be before it stops finding new
alternative splicing?

It does not stop, within the range tested. What falls is the rate, and where it
falls below a threshold depends on how strongly the gene is expressed. In the
heart cohort, the share of a category gained per extra 50M reads drops below 1%
between 100M and 150M for genes with 0.5 <= TPM < 10, at 150-200M for
0.1 <= TPM < 0.5, and not until 250-300M for the least expressed. For genes
above 10 TPM it never drops below 1% within 300M, but that category holds only
1307 genes against 13436 in the lowest, so the same handful of new detections
reads as a larger share; its absolute gain per step falls from 76 genes to 16.
Meanwhile each additional gene costs between one and a half and two and a half
times more than the one before, for every 50M step past 100M.

Where to draw the line is a trade between those two, per expression level.
`results/Table_3.tsv` and `results/Figure_7_values.tsv` hold the numbers;
`02_downsampled_cohorts` and `05_sequencing_cost` compute them. Percentages are
shares of `data/expression_category_totals.tsv`, whose denominators differ by an
order of magnitude between categories, so read the absolute columns alongside
them.

**Sequencing depth** always means mapped reads, never per-junction read support.
**Genes with AS** and **junctions involved in AS** are MAJIQ local splicing
variants with 0.05 < PSI < 0.95. The adipose, heart and hypothalamus samples
were each downsampled 100 times at every depth from 50M to 300M; a junction
found in all 100 subsamples is **robust**, in fewer **sporadic**. The
SARS-CoV-2 samples were not downsampled: each was sequenced once, at 60-150M.
An **expression category** is a TPM bin.

## Layout

```
notebooks/         the analysis, one notebook per question
asdepth/           readers for the MAJIQ summary tables, plot styling
scripts/           the steps too slow or too large to sit in a notebook
data/annotation/   per-gene type, length, compactness and GC content
data/derived/      the small tables the notebooks read
results/           figures and tables written by the notebooks
tests/             invariants the analysis has to satisfy
webapp/            source of sequencing-depth4splicing.com
```

Flow is one way: notebooks read `data/` and write `results/`, and nothing under
`results/` is read back in anywhere. Only file reading and plot styling live in
`asdepth`; the statistics stay in the notebooks.

Each notebook locates the top of the clone itself, so it runs from anywhere
inside it:

```
python3 -m nbconvert --to notebook --execute --inplace notebooks/04_gtex_and_tcga.ipynb
```

## The notebooks

| Notebook | Question | Download | Runs in |
|---|---|---|---|
| `01_sars_cov_2_cohort` | Detections per depth across samples sequenced 60-150M | yes | ~5 min |
| `02_downsampled_cohorts` | The same, downsampled 100x per depth to 300M | yes | ~15 min |
| `03_junction_features` | Which kinds of junction do the extra reads add? | no | seconds |
| `04_gtex_and_tcga` | How much do GTEx and TCGA miss at their usual depth? | no | seconds |
| `05_sequencing_cost` | What does each additional detection cost? | no | seconds |
| `06_detectability_model` | Which gene properties predict detection? | no | seconds |
| `07_biological_heterogeneity` | Does heterogeneity rather than depth drive discovery? | yes | ~1 min |
| `08_pathway_enrichment` | Are the pathways found only at 200M coherent? | no | seconds |

## Data

Five of the eight notebooks run straight after a clone. The other three read the
`*_summary.tsv` splicing tables, which are on Zenodo,
[10.5281/zenodo.11655945](https://doi.org/10.5281/zenodo.11655945), together
with the raw MAJIQ output; `07_biological_heterogeneity` also needs the
rnaseqmut output. Each notebook names the paths it expects in its first cell.

Two larger intermediates are not distributed. The junction feature table behind
`03_junction_features` and `06_detectability_model` is a gigabyte;
`scripts/build_junction_feature_table.py` builds it from the summary tables, a
GENCODE annotation and `data/annotation/`, and
`scripts/summarise_junction_features.py` reduces it to the four small tables in
`data/derived/` that those notebooks read. The DICAST-unified event tables that
NEASE was run on are not distributed either; its output is in
`data/derived/nease/`, which is what `08_pathway_enrichment` reads.

`data/README.md` says where every shipped file comes from.

## Building the inputs

```
bash    scripts/har_genes.sh                     # data/annotation/har_genes.txt
python3 scripts/build_junction_feature_table.py  # the junction feature table
python3 scripts/summarise_junction_features.py   # data/derived/feature_*.tsv
python3 scripts/fit_gee_models.py                # data/derived/gee_coefficients.tsv
```

`tpm.py` normalises featureCounts output to TPM.
`count_genes_from_majiq.py` and `count_junctions_from_majiq.py` each print one
line per MAJIQ output file, which is the format of
`data/derived/extrapolation_counts/`.
`make_depth_estimator_coeffs.py` builds `webapp/data/seq_depth_coeffs.json`.

## Checks

```
python3 -m pytest tests
```

Twelve invariants rather than stored numbers: the sub-categories of a feature
partition its junctions at every depth, the three expression tertiles have
slopes summing to zero, an FDR is never below its raw p-value, deeper
sequencing puts every cohort above its own baseline, junctions outnumber genes
and are cheaper to find, each additional detection costs more than the last, the
detectability model and the composition analysis agree on expression, and every
reported pathway is significant at 200M and not before. One check reads a
miniature summary file counted by hand. CI runs them on every push together with
the five notebooks that need no download.

## Requirements

Python 3.10 to 3.12 and `requirements.txt`; the pinned numpy does not build on
3.13.

MIT licensed, see [LICENSE](LICENSE). Citation metadata in `CITATION.cff`.
