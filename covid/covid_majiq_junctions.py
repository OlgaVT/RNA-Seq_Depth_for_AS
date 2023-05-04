import glob
from collections import Counter
import sys
import os

majiq_path = sys.argv[1] #a majiq tsv file
tpm_path = sys.argv[2] # a _TPM file
summary_path = sys.argv[3] # a STAR.Log.final.out file
t1 = float(sys.argv[4]) # a threshold for a low expression
t2 = float(sys.argv[5]) #a treshold for a high expression

exp_level = {}

tpm_file = open(tpm_path)
for line in tpm_file:
    try:
        if float(line.split()[-1]) >= t2:
            exp_level[line.split()[0]] = 'high'
        elif float(line.split()[-1]) < t1:
            exp_level[line.split()[0]] = 'low'
        else:
            exp_level[line.split()[0]] = 'mid'
    except:
        pass
tpm_file.close() #divide all genes according their expression levels

gene_set_low = set()
gene_set_mid = set()
gene_set_high = set()

majiq_file = open(majiq_path)
for line in majiq_file:
    if line.split()[0] != 'Gene':
        gene_id = line.split()[0].split(':')[1]
        junctions = line.split()[10].split(';')
        #gene_id = line.split()[0] #for the simulated data set
        if exp_level[gene_id] == 'low':
            for j in junctions:
                gene_set_low.add(gene_id + '_' + j)
        elif exp_level[gene_id] == 'mid':
            for j in junctions:
                gene_set_mid.add(gene_id + '_' + j)
        if  exp_level[gene_id] == 'high':
            for j in junctions:
                gene_set_high.add(gene_id + '_' + j)
majiq_file.close()

summary_file = open(summary_path)
for line in summary_file:
    if line.startswith('                          Number of input reads |'):
        seq_depth_1 = int(line.split()[-1])
    if line.startswith('                   Uniquely mapped reads number |'):
        seq_depth_2 = int(line.split()[-1])
summary_file.close()


print('high', len(gene_set_high), seq_depth_1, seq_depth_2, tpm_path)
print('mid', len(gene_set_mid), seq_depth_1, seq_depth_2, tpm_path)
print('low', len(gene_set_low), seq_depth_1, seq_depth_2, tpm_path)
