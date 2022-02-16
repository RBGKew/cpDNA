#!/bin/bash

#SBATCH --job-name="cpDNA"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=8000
ncpu=1

## Kew
# module load python/3.7.9
# module load amas/1.0
# module load iqtree/2.0.6
# module load fasttree
# module load trimal/1.4.1
# module load mafft/7.471
## JHI
source activate phylo

### Input ###
Label_file=$1 #cpTree0.1_labels.txt

### Parameters ### 
min_len_r=0
min_sample_r=0
slurmThrottle=80

### Derived input ###
project="ftp"
echo $project

### Concat all samples ### 
## Add all genes from all samples in a single file
mkdir -p $project
mkdir -p $project/Genes; mkdir -p $project/Aligned_genes; mkdir -p $project/Aligned_genes_trimmed
mkdir -p $project/by_sample_contigs; mkdir -p $project/by_sample_genes;
cd $project
> AllSamples_Allgenes.fasta
while read iline; do
	sample="$(cut -d' ' -f1 <<<"$iline")"
	label="$(cut -d' ' -f2 <<<"$iline")"
	cp ../CPgenes_by_sample/"$sample"_CPgenes.fasta by_sample_genes/"$sample"-genes.fasta
	cp ../Data/"$sample"_pt.fasta by_sample_contigs/"$sample"-contigs.fasta
	cat ../CPgenes_by_sample/"$sample"_CPgenes.fasta >> AllSamples_Allgenes.fasta
done < ../$Label_file

python ../scripts/Transpose_s2g.py -f AllSamples_Allgenes.fasta --min_len_r $min_len_r --min_sample_r $min_sample_r --out_dir Genes/ --genes_kept_file ls_genes.txt

### Align and trim genes ### 
Ngenes=$(wc -l ls_genes.txt | awk '{ print $1 }')
if (( $Ngenes > 0 )); then
	sbatch --wait -p medium --array=1-${Ngenes}%$slurmThrottle ../scripts/Align_gene.sh ls_genes.txt
fi




# http://www.iqtree.org/doc/Advanced-Tutorial
# iqtree -s concat.fasta -p concat.partitions -m GTR+F+I+G4 --prefix test_concat -T 4
#iqtree -p Aligned_genes_trimmed -m MFP --prefix test_p -B 1000 -T AUTO
# iqtree -s Concat_aligned_38genes.fasta -m MFP -B 1000 -T AUTO