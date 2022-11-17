import glob
from collections import Counter
import sys
import os

majiq_path = sys.argv[1] #path to a majiq tsv file
tpm_path = sys.argv[2] #path to a TPM file
summary_path = sys.argv[3] #path to a STAR Log.final.out file
t1 = float(sys.argv[4]) # a threshold for a low expression
t2 = float(sys.argv[5]) # a treshold for a high expression

exp_level = {}

tpm_file = open(tpm_path)
for line in tpm_file:
    try:
        if float(line.split()[-2]) >= t2:
            exp_level[line.split()[0]] = 'high'
        elif float(line.split()[-2]) < t1:
            exp_level[line.split()[0]] = 'low'
        else:
            exp_level[line.split()[0]] = 'mid'
    except:
        pass
tpm_file.close() #divide all genes according their expression levels

gene_list_low = set()
gene_list_mid = set()
gene_list_high = set()

majiq_file = open(majiq_path)
for line in majiq_file:
    if line.split()[0] != 'Gene':
        gene_id = line.split()[0].split(':')[1]
        #gene_id = line.split()[0] #for the simulated data set
        if exp_level[gene_id] == 'low':
            gene_list_low.add(gene_id)
        elif exp_level[gene_id] == 'mid':
            gene_list_mid.add(gene_id)
        if  exp_level[gene_id] == 'high':
            gene_list_high.add(gene_id)
majiq_file.close()

summary_file = open(summary_path)
for line in summary_file:
    if line.startswith('                          Number of input reads |'):
        seq_depth_1 = int(line.split()[-1])
    if line.startswith('                   Uniquely mapped reads number |'):
        seq_depth_2 = int(line.split()[-1])
summary_file.close()


print('high', len(gene_list_high), seq_depth_1, seq_depth_2, tpm_path)
print('mid', len(gene_list_mid), seq_depth_1, seq_depth_2, tpm_path)
print('low', len(gene_list_low), seq_depth_1, seq_depth_2, tpm_path)
