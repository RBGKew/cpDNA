#!/usr/bin/env python
# coding: utf-8

# In[120]:


import pandas as pd
from Bio import SeqIO
import argparse


# In[121]:


# params = {}
# params['filename'] = '../AllSamples_Allgenes.fasta'
# params['min_len_r'] = 0; params['min_sample_r'] = 0; 
# params['out_dir'] = '../Genes/'
# params['genes_kept_file'] = 'ls_genes_kept.txt'
# params['acc_genes'] = True
# params['rename'] = '../cpTree_v7/cpTree_v7_tree_DataDistribution.txt'
# params['rename'] = False


# In[122]:


parser = argparse.ArgumentParser(description='Get genes from all samples')
parser.add_argument("-f", help="fasta file", dest="filename", required=True)
parser.add_argument("--min_len_r", help="min_len_r", action="store", type=float, default=0.2)
parser.add_argument("--min_sample_r", help="min_sample_r", action="store", type=float, default=0.5)
parser.add_argument("--acc_genes", help="Keep only accepted genes (defined in script)", action="store_true", default=False)
parser.add_argument("--out_dir", help="Output dir", action="store", default='Genes/')
parser.add_argument("--rename", help="Rename Samples. 2nd column is new name, 3rd is label", action="store", default=False)
parser.add_argument("--genes_kept_file", help="Output genes kept", action="store", default='ls_genes.txt')
params = vars(parser.parse_args());


# In[123]:


print('Parameters:',params)


# In[124]:


# records
print('Loading records...')
rec_dc = {}
Records = SeqIO.to_dict(SeqIO.parse(params['filename'],format='fasta'))
for id, record in Records.items():
    rec_dc[id] = len(record.seq)


# In[125]:


rec_df = pd.DataFrame.from_dict(rec_dc,orient='index').reset_index().rename(columns={'index':'seqid',0:'len'})
rec_df[['Sample_Name','gene']] = rec_df.seqid.str.split('-',expand=True)
print('found',rec_df.shape[0],'sequences,',rec_df.Sample_Name.nunique(),'samples,',rec_df.gene.nunique(),'genes')


# In[126]:


if params['rename']!=False:
    rename_df = pd.read_table(params['rename'],sep=' ',header=None)
    rename_df.columns = ['Sample_Name','Public_Name','Label']
    rec_df = pd.merge(rec_df,rename_df,how='left',on='Sample_Name')
    print(rec_df.isna().sum().to_dict())


# In[127]:


rec_df.to_csv(params['filename'].replace('.fasta','_All_SeqTable.csv'),index=False)


# In[128]:


if params['acc_genes']==True:
    accepted_genes = ['accD','atpA','atpB','atpE','atpF','atpH','atpI','ccsA','cemA','clpP','matK','ndhA','ndhB','ndhC','ndhD','ndhE',
                  'ndhF','ndhG','ndhH','ndhI','ndhJ','ndhK','petA','petB','petD','petG','petL','petN','psaA','psaB','psaC','psaI',
                  'psaJ','psbA','psbB','psbC','psbD','psbE','psbF','psbH','psbI','psbJ','psbK','psbL','psbM','psbN','psbT','psbZ',
                  'rbcL','rpl2','rpl14','rpl16','rpl20','rpl22','rpl23','rpl32','rpl33','rpl36','rpoA','rpoB','rpoC1','rpoC2','rps2',
                  'rps3','rps4','rps7','rps8','rps11','rps12','rps14','rps15','rps16','rps18','rps19','ycf3','ycf4']
    print('keeping only',len(accepted_genes),'accepted genes')
    rec_df = rec_df[rec_df.gene.isin(accepted_genes)]
    print('>',rec_df.shape[0],'sequences,',rec_df.Sample_Name.nunique(),'samples,',rec_df.gene.nunique(),'genes')


# In[129]:


# Median Sum of gene length for good recoveries
Stats_len = rec_df.groupby('Sample_Name').agg({'len':'sum','gene':'count'})#.to_frame()
print('Median gene count:',Stats_len.gene.median(),'Median sum of gene lengths:',Stats_len.len.median())
Stats_len_good = Stats_len[Stats_len.gene>=rec_df.gene.nunique()*0.9]
median_sumlen = Stats_len_good.len.median()
print(round(Stats_len_good.shape[0]/Stats_len.shape[0]*100,0),'% of samples have',rec_df.gene.nunique(),
      'genes. Median Sum of genes length:',median_sumlen)


# In[130]:


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


# In[131]:


# Update gene list
genes_df = rec_df.groupby('gene').size().to_frame().reset_index()
genes_df.gene.to_csv(params['genes_kept_file'],index=False,header=None)


# In[132]:


# Write output, gene by gene
for idx, row in genes_df.iterrows():
    rec_gene = []
    rec_gene_df = rec_df[rec_df.gene==row.gene]
    for idx, row in rec_gene_df.iterrows():
        if params['rename']!=False:
            Records[row.seqid].id = row.Public_Name
            Records[row.seqid].description = row.Label + ' ' + Records[row.seqid].description
        else:
            Records[row.seqid].id = row.Sample_Name
        Records[row.seqid].description = Records[row.seqid].description.replace('reference_coverage','reference_overlap').replace('organism-gene','Reference_organism-gene')
        rec_gene.append(Records[row.seqid])
    SeqIO.write(rec_gene, params['out_dir'] + row.gene + '.fasta',format='fasta')

