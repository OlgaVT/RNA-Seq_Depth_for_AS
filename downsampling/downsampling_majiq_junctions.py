import glob
from collections import Counter
import sys
import os

majiq_path = sys.argv[1]  #a PATH to majiq tsv files with majiq results per sample in a glob format, e.g. "/path/to/majiq_tsv/*"
tpm_path = sys.argv[2]  # a TPM file
t1 = float(sys.argv[3])  # a low expression threshold 
t2 = float(sys.argv[4])  # a high expression threshold
output = sys.argv[5]  #  "_x_junecions_in_y_samples"

if os.path.exists(output):
    os.remove(output)

majiq_file_list = glob.glob(majiq_path)

exp_level = {}

tpm_file = open(tpm_path)
for line in tpm_file:
    try:
        if float(line.split()[-1]) >= t2:
            exp_level[line.split()[1]] = 'high'
        elif float(line.split()[-1]) < t1:
            exp_level[line.split()[1]] = 'low'
        else:
            exp_level[line.split()[1]] = 'mid'
    except:
        pass
tpm_file.close() #divide all gene according their expression levels

gene_list_low = []
gene_list_mid = []
gene_list_high = []


for i in majiq_file_list:
    majiq_file = open(i)
    print(i)
    gene_set_high = set()
    gene_set_low = set()
    gene_set_mid = set()
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
    gene_list_low += list(gene_set_low)
    gene_list_mid += list(gene_set_mid)
    gene_list_high += list(gene_set_high)

def Counting(gene_list, explevel):
    genes_counter = Counter()
    number_genes_counter = Counter()
    for i in gene_list:
        genes_counter[i] += 1

    for k,v in genes_counter.items():
        number_genes_counter[v] += 1

    for k,v in number_genes_counter.items():
        x_junctions_in_y_samples = open(sys.argv[5], 'a')
        result = '\t'.join((str(k), str(v), explevel))
        x_junctions_in_y_samples.write(result + '\n')
        x_junctions_in_y_samples.close()

Counting(gene_list_low, 'low')
Counting(gene_list_mid, 'mid')
Counting(gene_list_high, 'high')
