#!/bin/bash

#SBATCH --job-name="cpDNA"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=20
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --mem=80000
ncpu=20

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
# project=$2

### Parameters ### 
Reference=References/cpDNA_RefCDS.fasta
min_cov=50
min_len_r=0.1
min_sample_r=0.5
slurmThrottle=80

### Derived input ###
project="$(echo $Label_file | sed 's/_labels.txt//g')"_"$min_len_r"
echo $project

### Extract cp genes from contigs ### 
mkdir -p Data_clean; mkdir -p Blastn; mkdir -p CPgenes_by_sample
## List remaining samples to do
samples_list=tmp_samples.txt
python scripts/list_blast_samples.py -i $Label_file -o $samples_list
## get cp genes
Nsamples=$(wc -l $samples_list | awk '{ print $1 }')
if (( $Nsamples > 0 )); then
	sbatch --wait -p short --array=1-${Nsamples}%$slurmThrottle scripts/Get_CPgenes.sh $Reference $samples_list $min_cov
fi

### Concat and filter all samples ### 
## Add all genes from all samples in a single file
mkdir -p $project
mkdir -p $project/Genes; mkdir -p $project/Aligned_genes; mkdir -p $project/Aligned_genes_trimmed
cd $project
> AllSamples_Allgenes.fasta
while read iline; do
	sample="$(cut -d' ' -f1 <<<"$iline")"
	cat ../CPgenes_by_sample/"$sample"_CPgenes.fasta >> AllSamples_Allgenes.fasta
done < ../$Label_file

python ../scripts/Transpose_s2g.py -f AllSamples_Allgenes.fasta --min_len_r $min_len_r --min_sample_r $min_sample_r --out_dir Genes/ --genes_kept_file ls_genes.txt --acc_genes

### Align and trim genes ### 
Ngenes=$(wc -l ls_genes.txt | awk '{ print $1 }')
if (( $Ngenes > 0 )); then
	sbatch --wait -p medium --array=1-${Ngenes}%$slurmThrottle ../scripts/Align_gene.sh ls_genes.txt
fi

### Concat genes ###
# https://github.com/marekborowiec/AMAS
# AMAS.py concat -f fasta -d dna -i Aligned_genes_trimmed/*.fasta -c 4 -u fasta --part-format raxml --concat-part concat.partitions --concat-out concat.fasta
# sed -i 's/?/-/g' concat.fasta # Replace missing genes (?) by insertions (-)

### Build tree ###
# fasttree -gtr -nt -gamma -log "$project"_fasttree.log < concat.fasta > "$project"_fasttree.nwk
# nw_rename -l "$project"_fasttree.nwk ../$Label_file > "$project"_fasttree_labelled.nwk
date +"%T"
iqtree -s Aligned_genes_trimmed/ -m GTR+G -B 1000 -T $ncpu --prefix "$project" --seqtype DNA -v
date +"%T"
nw_rename -l "$project".treefile ../$Label_file > "$project"_iqtree_labelled.nwk



# http://www.iqtree.org/doc/Advanced-Tutorial
# iqtree -s concat.fasta -p concat.partitions -m GTR+F+I+G4 --prefix test_concat -T 4
#iqtree -p Aligned_genes_trimmed -m MFP --prefix test_p -B 1000 -T AUTO
# iqtree -s Concat_aligned_38genes.fasta -m MFP -B 1000 -T AUTO