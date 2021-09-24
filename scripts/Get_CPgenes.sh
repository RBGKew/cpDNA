#!/bin/bash
#SBATCH --job-name="Get_CPgenes"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=8000
ncpu=4

#JHI
source activate phylo

Reference=$1
Samples_file=$2
min_cov=$3

iline=$(sed -n "$SLURM_ARRAY_TASK_ID"p $Samples_file)
echo $iline
sample="$(cut -d',' -f1 <<<"$iline")"

python scripts/copy_fasta_cutid.py Data/"$sample"_pt.fasta Data_clean/"$sample"_pt.fasta
makeblastdb -in Data_clean/"$sample"_pt.fasta -dbtype nucl -parse_seqids

blastn -query $Reference -db Data_clean/"$sample"_pt.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore sseq" -out Blastn/"$sample"_CPgenes.blastn -max_target_seqs 5 -num_threads $ncpu

python scripts/targets_from_blast.py Blastn/"$sample"_CPgenes.blastn CPgenes_by_sample/"$sample"_CPgenes.fasta $min_cov

