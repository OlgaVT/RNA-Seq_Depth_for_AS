import glob
from collections import Counter
import sys
import os

majiq_path = sys.argv[1] #path to a majiq tsv file
seq_depth =  sys.argv[2]
seq_depth_unique =  sys.argv[3]

gene_set = set()
majiq_file = open(majiq_path)
for line in majiq_file:
    if line.split()[0] != 'Gene':
        gene_id = line.split()[0].split('.')[0]
        gene_set.add(gene_id)
majiq_file.close()


print(majiq_path, len(gene_set), seq_depth, seq_depth_unique)
