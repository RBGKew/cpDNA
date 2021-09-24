from Bio import SeqIO
import sys
in_fasta = sys.argv[1]
out_fasta = sys.argv[2]

max_lenid=20
records = list(SeqIO.parse(in_fasta, "fasta"))
for i in range(len(records)):
	if len(records[i].id)>max_lenid:
		records[i].id = records[i].id[:max_lenid].replace(',','_')

SeqIO.write(records,out_fasta,format='fasta')

# python scripts/copy_fasta_cutid.py ../GetOrganelles/PAFTOL/fasta_pt/PAFTOL_001191_pt.fasta fasta_pt/PAFTOL_001191_pt.fasta