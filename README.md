# RNA-Seq_Depth_for_AS

This is a repository for the paper 'RNA sequencing depth guidelines for the study of alternative splicing' that contains all codes for visualization and the part of the processed data. The rest of the data: outputs of the MAJIQ analysis and summary files could be downloaded here:

Input files: Zenodo [10.5281/zenodo.11655945]

1) Figure1: the visualization Jupyter notebook for Figure 1 in the part 'SARS-CoV-2 cohort'. The covid_summary.tsv is at the Zenodo.
The columns of the covid_summary.tsv:
Phenotype Sample Gene Junction PSI Seq_Depth TPM 

3) Figure3: the visualization Jupyter notebook for Figure 3 in the part 'Deep-sequenced cohorts'. The X_summary.tsv files are at the Zenodo.
The columns of the X_summary.tsv:
Phenotype Sample Gene Junction PSI Seq_Depth TPM

5) Enrichment (For the section "The biological relevance of AS events found at 200M reads"): the visualization Jupyter notebook for NEASE enrichment analysis. The /Enrichment/Files folder contains the following files:
   - X_more_than_X.tsv - genes detected exclusively at 200M reads
   - _unified - the coordinates of the detected events unified using the DICAST
   - _nease - the results of the NEASE enrichment analysis
6) Figure4: the visualization Jupyter notebook and the input files for Figure 4 in the part 'Comparison of GTEx and TCGA cohorts with the deep-sequenced cohorts'.
7) Figure 5: the visualization Jupyter notebook and the input files for Figure 4 in the part 'The cost of AS detection in the deep-sequenced cohorts'.
8) TableS7: the Jupyter notebook and the input files to calculate the amount of novel junctions
9) Others: TPM.py script for normalizing the raw counts from featureCounts output.
