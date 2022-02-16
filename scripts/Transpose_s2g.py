#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from Bio import SeqIO
import argparse


# In[7]:


# params = {}
# params['filename'] = '../allgenes_allsamples.fasta'
# params['min_len_r'] = 0.2; params['min_sample_r'] = 0.5; 
# params['out_dir'] = '../Genes/'
# params['genes_kept_file'] = 'ls_genes_kept.txt'
# params['acc_genes'] = True


# In[4]:


parser = argparse.ArgumentParser(description='Get genes from all samples')
parser.add_argument("-f", help="fasta file", dest="filename", required=True)
parser.add_argument("--min_len_r", help="min_len_r", action="store", type=float, default=0.2)
parser.add_argument("--min_sample_r", help="min_sample_r", action="store", type=float, default=0.5)
parser.add_argument("--acc_genes", help="Keep only accepted genes (defined in script)", action="store_true", default=False)
parser.add_argument("--out_dir", help="Output dir", action="store", default='Genes/')
parser.add_argument("--genes_kept_file", help="Output genes kept", action="store", default='ls_genes.txt')
params = vars(parser.parse_args());


# In[3]:


print('Parameters:',params)


# In[4]:


# records
print('Loading records...')
rec_dc = {}
Records = SeqIO.to_dict(SeqIO.parse(params['filename'],format='fasta'))
for id, record in Records.items():
    rec_dc[id] = len(record.seq)


# In[5]:


rec_df = pd.DataFrame.from_dict(rec_dc,orient='index').reset_index().rename(columns={'index':'seqid',0:'len'})
rec_df[['Sample_Name','gene']] = rec_df.seqid.str.split('-',expand=True)
rec_df.to_csv(params['filename'].replace('.fasta','_All_SeqTable.csv'),index=False)
print('found',rec_df.shape[0],'sequences,',rec_df.Sample_Name.nunique(),'samples,',rec_df.gene.nunique(),'genes')


# In[6]:


if params['acc_genes']==True:
    accepted_genes = ['accD','atpA','atpB','atpE','atpF','atpH','atpI','ccsA','cemA','clpP','matK','ndhA','ndhB','ndhC','ndhD','ndhE',
                  'ndhF','ndhG','ndhH','ndhI','ndhJ','ndhK','petA','petB','petD','petG','petL','petN','psaA','psaB','psaC','psaI',
                  'psaJ','psbA','psbB','psbC','psbD','psbE','psbF','psbH','psbI','psbJ','psbK','psbL','psbM','psbN','psbT','psbZ',
                  'rbcL','rpl2','rpl14','rpl16','rpl20','rpl22','rpl23','rpl32','rpl33','rpl36','rpoA','rpoB','rpoC1','rpoC2','rps2',
                  'rps3','rps4','rps7','rps8','rps11','rps12','rps14','rps15','rps16','rps18','rps19','ycf3','ycf4']
    print('keeping only',len(accepted_genes),'accepted genes')
    rec_df = rec_df[rec_df.gene.isin(accepted_genes)]
    print('>',rec_df.shape[0],'sequences,',rec_df.Sample_Name.nunique(),'samples,',rec_df.gene.nunique(),'genes')


# In[17]:


# Median Sum of gene length for good recoveries
Stats_len = rec_df.groupby('Sample_Name').agg({'len':'sum','gene':'count'})#.to_frame()
Stats_len_good = Stats_len[Stats_len.gene>=rec_df.gene.nunique()]
median_sumlen = Stats_len_good.len.median()
print(round(Stats_len_good.shape[0]/Stats_len.shape[0]*100,0),'% of samples have',rec_df.gene.nunique(),
      'genes. Median Sum of genes length:',median_sumlen)


# In[24]:


# Filter by sum of contigs length
min_contigs_len = median_sumlen*params['min_len_r']
print('min_contigs_len:',min_contigs_len, (Stats_len['len']<min_contigs_len).sum(),'samples failed')
samples_pass = list(Stats_len[Stats_len['len']>=min_contigs_len].index)
rec_df = rec_df[rec_df.Sample_Name.isin(samples_pass)]

# Filter by number of samples per gene
min_sample = round(params['min_sample_r']*rec_df.Sample_Name.nunique())
tmp = rec_df.groupby('gene').size().to_frame().rename(columns={0:'N'})
print('min_sample:',min_sample, (tmp.N<min_sample).sum(),'genes failed')
genes_pass = list(tmp[tmp.N>=min_sample].index)
rec_df = rec_df[rec_df.gene.isin(genes_pass)]

print('after filtering',rec_df.shape[0],'sequences,',rec_df.Sample_Name.nunique(),'samples,',rec_df.gene.nunique(),'genes')
rec_df.to_csv(params['filename'].replace('.fasta','_Sel_SeqTable.csv'),index=False)


# In[26]:


# Update gene list
genes_df = rec_df.groupby('gene').size().to_frame().reset_index()
genes_df.gene.to_csv(params['genes_kept_file'],index=False,header=None)


# In[28]:


# Write output, gene by gene
for idx, row in genes_df.iterrows():
    rec_gene = []
    seqid_ls = list(rec_df[rec_df.gene==row.gene].seqid)
    for seqid in seqid_ls:
        Records[seqid].id = Records[seqid].id.split('-')[0]
        rec_gene.append(Records[seqid])
    SeqIO.write(rec_gene, params['out_dir'] + row.gene + '.fasta',format='fasta')

