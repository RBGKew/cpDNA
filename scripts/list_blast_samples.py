#!/usr/bin/env python
# coding: utf-8

# Example: python scripts/list_blast_samples.py -i R1_labels.txt -o tmp_samples.txt


import pandas as pd
import argparse; import os


parser = argparse.ArgumentParser()
parser.add_argument("-i"); parser.add_argument("-o")
opts = parser.parse_args()
label_file = opts.i;  out_file = opts.o

db = pd.read_table(label_file,sep=' ',header=None).rename(columns={0:'Sample',1:'Label'})

# List existing cpGenes files cards output list of samples to blast
samples_done = [filename.replace('_CPgenes.fasta','') for filename in os.listdir('CPgenes_by_sample/')]
print('\nfound',len(samples_done),'sample gene fasta')

# Output list of samples to blast
samples_todo = db[db.Sample.isin(samples_done)==False]
print(samples_todo.shape[0],'samples to blast, ',db[db.Sample.isin(samples_done)].shape[0],'samples done')
samples_todo.Sample.to_csv(out_file,index=False,header=False)

