#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#python scripts/targets_from_blast.py Blastn/PAFTOL_012529_cpRefGenes.blastn cpGenes_by_sample/PAFTOL_012529_cpGenes.fasta 50


# In[1]:


import sys
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq; from Bio.SeqRecord import SeqRecord


# In[2]:


def blast_results(blast_file):
    #load blast result
    sample_blast=pd.read_table(blast_file, sep='\t', header=None) 
    
    # Because references are blasted against the contigs, I'm swapping query and subject fields to make the analysis clearer.
    qseqid='sseqid'; sseqid='qseqid'; qlen='slen'; qstart='sstart'; 
    qend='send'; slen='qlen'; sstart='qstart'; send='qend'
    sample_blast.columns=[qseqid, sseqid, 'pident', 'length', 'mismatch', 'gapopen', qlen, qstart, 
                              qend, slen, sstart, send, 'evalue', 'bitscore','qseq']
    # Now query is the contigs and subject is the references

    #split reference id to get ref species and ref gene
    sample_blast[['ref_sp', 'starget']] = sample_blast['sseqid'].str.split('-', n=1, expand=True)
    
    # Get direction of match and reverse complement if necessary
    sample_blast['qdir'] = (sample_blast.qstart<sample_blast.qend).replace({True:'For',False:'Rev'})
    sample_blast['sdir'] = (sample_blast.sstart<sample_blast.send).replace({True:'For',False:'Rev'})
    for idx, row in sample_blast.iterrows():
        if (row.qdir=='Rev'):
            sample_blast.loc[idx,'qseq'] = str(Seq(row['qseq']).reverse_complement())
            tmp = row['qstart']; sample_blast.loc[idx,'qstart'] = row['qend']; sample_blast.loc[idx,'qend']=tmp
        elif (row.sdir=='Rev'):
            sample_blast.loc[idx,'qseq'] = str(Seq(row['qseq']).reverse_complement())
            tmp = row['sstart']; sample_blast.loc[idx,'sstart'] = row['send']; sample_blast.loc[idx,'send']=tmp
    
    #calculate overlap
    sample_blast['scover']=round((sample_blast.send-sample_blast.sstart+1)/sample_blast.slen*100,1)
    sample_blast['qcover']=round((sample_blast.qend-sample_blast.qstart+1)/sample_blast.qlen*100,1)
    sample_blast['overlap']=False
    
    return sample_blast[['qseqid', 'sseqid', 'starget', 'ref_sp', 'pident', 'evalue', 
                         'bitscore','length', 'qlen', 'slen', 'qstart', 'qend', 'sstart','send','qdir','sdir',
                        'scover', 'qcover', 'overlap', 'qseq']]


# In[3]:


# Make dictionary of nucleotides
def seq_str2dc(seq_dc):
    return dict(zip(list(np.array([i for i in range(seq_dc['sstart'],seq_dc['send']+1)])), seq_dc['qseq']))

def overlap_resolve(df):
    df = df.sort_values('bitscore',ascending=False)
    best = df.to_dict('records')[0] # Take best bitscore as backbone
    best['seq_dc'] = seq_str2dc(best)
    current_length = best['length'];
    # Iterate over other matches to complete the reference (s)
    for idx, row in df.iterrows():
        row_dc = row.to_dict()
        row_dc['seq_dc'] = seq_str2dc(row_dc)
        best['seq_dc'] = {**row_dc['seq_dc'], **best['seq_dc']} # Add only new positions
        # update bitscore if length as increased
        best['length'] = len(best['seq_dc'])
        if best['length']>current_length:
            update_len =  best['length'] - current_length
            update_ratio = update_len/row_dc['length']
            best['bitscore'] = round(best['bitscore'] + (update_ratio*row_dc['bitscore']),0)
            best['pident'] = round( ((best['pident']*(current_length)) + (update_len*row_dc['pident'])) / best['length'],0)
            best['scover'] = round((best['length'])/best['slen']*100,1)
            best['overlap']=True
            current_length = best['length'];
    best['qseq'] = "".join([N[1] for N in sorted(best['seq_dc'].items())])
    del best["seq_dc"]
    return best

# test_df = pd.DataFrame([{'Sample':'seq3','pident':90,'bitscore':10,'sstart':1,'send':6,'slen':25,'qseq':'AGAGAG'},
#                         {'Sample':'seq1','pident':80,'bitscore':50,'sstart':10,'send':20,'slen':25,'qseq':'ACACACACACA'},
#                        {'Sample':'seq4','pident':80,'bitscore':9,'sstart':19,'send':25,'slen':25,'qseq':'XXXXXXX'},
#                         {'Sample':'seq2','pident':90,'bitscore':25,'sstart':20,'send':25,'slen':25,'qseq':'CTCTCT'}])
# test_df['length'] = test_df.qseq.str.len()
# test_df['scover'] = test_df.length/test_df.slen*100
# best = overlap_resolve(test_df)
# print(test_df)
# print(best)


# In[4]:


def writefasta(df, fasta_out, sample_name):
    out_seqs = []
    for index, iseq in df.iterrows():
        id = sample_name + '-' + iseq.starget
        description = '{organism-gene:' + iseq.ref_sp + '-' + iseq.starget + ', originalID:'         + iseq.qseqid.replace(',','') + ', length:' + str(iseq.length) + ', pident:' + str(round(iseq.pident,2)) + ', bitscore:'         + str(iseq.bitscore) + ', reference_coverage:' + str(iseq.scover) + ', overlap:' + str(iseq.overlap) + '}'
#         + ', sstart:' + str(iseq.sstart) + ', send:' + str(iseq.send) + ', qstart:' + str(iseq.qstart) + ', qend:' + str(iseq.qend) + '}'
        # + ';cut2codon:' + str(iseq.cut2codon) 
        
        # Remove insertions
        seq = Seq(iseq.qseq.replace('-', ''))
        record = SeqRecord(seq=seq,id=id,description=description)
        
        out_seqs.append(record)
        
    # Write sequences
    SeqIO.write(out_seqs,fasta_out,format='fasta')


# In[6]:


# Main
if __name__ == "__main__":
    input_blast_file=sys.argv[1]
    fasta_out=sys.argv[2]
#     input_blast_file= '../tmp_cpDNA/' + 'PAFTOL_006071_CPgenes.blastn'
#     fasta_out= '../tmp_cpDNA/' + 'PAFTOL_006071_CPgenes.fasta'
    
    sample_name=input_blast_file.split('/')[-1].replace('_CPgenes.blastn','')
    min_cov=50
    min_pid=70
    resolve_overlap=True
    restrict_genes=True
    
    
    print('\n',sample_name, input_blast_file, fasta_out)
    
    if os.path.exists(input_blast_file):
        if os.stat(input_blast_file).st_size != 0:
            # Load blast output
            df_blast=blast_results(blast_file=input_blast_file)
            print('All entries:',df_blast.shape[0], df_blast.starget.nunique())
            print(df_blast.groupby(['qdir','sdir']).size())
            
            ### Resolve overlaps: genes with multiple exons AND (if restrict_genes==False) genes with multiple matching contigs
            if resolve_overlap:
                genes_ex = ['atpF', 'clpP1', 'ndhA', 'ndhB', 'pafI', 'petB', 'petD', 'rpl16',
                   'rpl2', 'rpoC1', 'rps12', 'rps16', 'trnC_ACA', 'trnG_UCC', 'trnK_UUU',
                   'trnL_UAA', 'trnS_CGA', 'trnV_UAC']

                # sseqid with multiple entries
                sseqid_ex = df_blast.groupby('sseqid').size().to_frame().rename(columns={0:'N'}).reset_index();
                sseqid_ex = sseqid_ex[sseqid_ex.N>1]
                if sseqid_ex.shape[0]>0:
                    sseqid_ex[['ref_sp', 'starget']] = sseqid_ex['sseqid'].str.split('-', n=1, expand=True)
                    print(sseqid_ex.shape[0],'sseqid with >1 match',sseqid_ex.starget.nunique(),'genes,',
                          sseqid_ex[sseqid_ex.starget.isin(genes_ex)].starget.nunique(),'known split genes')
                    if restrict_genes:
                        sseqid_ex = sseqid_ex[sseqid_ex.starget.isin(genes_ex)]
                        
                # Resolve overlap
                if sseqid_ex.shape[0]>0:
                    print('restricted:',sseqid_ex.shape[0],'sseqid with >1 match in',sseqid_ex.starget.nunique(),'genes')
                    ls_overlap = []
                    for sseqid in sseqid_ex.sseqid:
                        ls_overlap.append(overlap_resolve(df_blast[df_blast.sseqid==sseqid]))
                    df_blast = pd.concat([df_blast, pd.DataFrame(ls_overlap)])
            
            ### Keep best match for each gene (starget)
            df_best = df_blast.sort_values('bitscore',ascending=False).groupby('starget').head(1).reset_index(drop=True)
            df_best = df_best[(df_best.scover>=min_cov) & (df_best.pident>=min_pid)]
            print('After filtering:',df_best.shape[0], df_best.starget.nunique(), df_best.length.sum(),
                 'overlap:',df_best.groupby('overlap').size().to_dict())
                                             
            if df_best.shape[0]>0:
                # Write fasta
                writefasta(df = df_best, fasta_out = fasta_out, sample_name = sample_name)

