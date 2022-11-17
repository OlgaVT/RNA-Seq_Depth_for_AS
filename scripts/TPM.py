import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep='\t', skiprows=1) #input from featurecounts
new_title = sys.argv[2]
#usage: python TPM.py featurecounts_output new_title

df.columns = ['Geneid', 'Chr', "Start", 'End', 'Strand', 'Length', new_title]
#TPM
df['RPK'] = df[new_title] * 1000 / df['Length']
per_million = df['RPK'].sum() / 1000000
df['TPM'] = df['RPK'] / per_million
#To file
df_result = df.drop(columns=['RPK'])
name = new_title + '_TPM'
df_result.to_csv(name, sep='\t')
