#!/bin/bash
#SBATCH --job-name="Get_CPgenes"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=8000
ncpu=1

genes_file=$1

gene=$(sed -n "$SLURM_ARRAY_TASK_ID"p $genes_file)
echo $gene

mafft --adjustdirection --auto Genes/"$gene".fasta > Aligned_genes/"$gene"_aln.fasta
sed -i 's/_R_//g' Aligned_genes/"$gene"_aln.fasta
trimal -in Aligned_genes/"$gene"_aln.fasta -out Aligned_genes_trimmed/"$gene"_alnT.fasta -gappyout


