#!/bin/bash

#SBATCH --job-name="cpDNA_pipeline"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=8000

### Commands
# ./cpDNA_pipeline.sh cpTree_v7_tree_DataDistribution.txt "distrib" 0
# ./cpDNA_pipeline.sh Chrysobalanaceae_cpTree_v7_clean_labels.txt "tree" 0.1

## Kew
# module load python/3.7.9
# module load amas/1.0
# module load iqtree/2.0.6
# module load fasttree
# module load trimal/1.4.1
# module load mafft/7.471
## JHI
###source activate phylo	# Paul B. - changed to set env outside of this script

### Input ###
Label_file=$1 #cpTree0.1_labels.txt
run_mode=$2 # filter, tree, distrib
min_len_r=$3

### Parameters ### 
Reference=References/cpDNA_RefCDS.fasta
min_cov=50
min_sample_r=0.5
slurmThrottle=80

### Derived input ###
if [ $run_mode = "distrib" ]; then 
	project="$(echo $Label_file | sed 's/.txt//g')"
else
	project="$(echo $Label_file | sed 's/_labels.txt//g')"_"$min_len_r"	
fi

echo $Label_file $project $run_mode $min_len_r

### Extract cp genes from contigs ###
if [ $run_mode = "filter" ]; then 
	"List samples for plastid gene recovery"
	mkdir -p Data_clean; mkdir -p Blastn; mkdir -p CPgenes_by_sample
	## List remaining samples to do
	samples_list=tmp_samples.txt
	python scripts/list_blast_samples.py -i $Label_file -o $samples_list
	## get cp genes
	Nsamples=$(wc -l $samples_list | awk '{ print $1 }')
	echo $Nsamples
	if (( $Nsamples > 0 )); then
		sbatch --wait -p short --array=1-${Nsamples}%$slurmThrottle scripts/Get_CPgenes.sh $Reference $samples_list $min_cov
	fi
fi


### Concat and filter all samples ### 
## Add all genes from all samples in a single file
mkdir -p $project
mkdir -p $project/Genes; 
cd $project
echo "Concat all gene sequences in one fasta file"
> AllSamples_Allgenes.fasta
while read iline; do
	sample="$(cut -d' ' -f1 <<<"$iline")"
	cat ../CPgenes_by_sample/"$sample"_CPgenes.fasta >> AllSamples_Allgenes.fasta
done < ../$Label_file

echo "Transpose samples to genes"
if [ $run_mode = "distrib" ]; then
	python ../scripts/Transpose_s2g.py -f AllSamples_Allgenes.fasta --min_len_r 0 --min_sample_r 0 --out_dir Genes/ --genes_kept_file ls_genes.txt --rename ../$Label_file
else
	python ../scripts/Transpose_s2g.py -f AllSamples_Allgenes.fasta --min_len_r $min_len_r --min_sample_r $min_sample_r --out_dir Genes/ --genes_kept_file ls_genes.txt --acc_genes
fi


### Align and trim genes ### 
if [ $run_mode = "tree" ]; then
	echo "Align and trim genes"
	mkdir -p Aligned_genes; mkdir -p Aligned_genes_trimmed
	Ngenes=$(wc -l ls_genes.txt | awk '{ print $1 }')
	if (( $Ngenes > 0 )); then
		sbatch --wait -p medium --array=1-${Ngenes}%$slurmThrottle ../scripts/Align_gene.sh ls_genes.txt
	fi
fi


### Build tree ###
if [ $run_mode = "tree" ]; then
	echo "Launch tree building job"
	sbatch -p long ../scripts/tree_building.sh $project $Label_file
fi

### Data Distribution (Dryad) ###
if [ $run_mode = "distrib" ]; then
	echo "Copying data in Data Distribution folder"
	cp ../$Label_file .
	### Copy contigs and genes by sample
	mkdir -p Contigs; mkdir -p Genes_By_Sample;
	while read iline; do
		Sample="$(cut -d' ' -f1 <<<"$iline")"
		PublicName="$(cut -d' ' -f2 <<<"$iline")"
		cp ../CPgenes_by_sample/"$Sample"_CPgenes.fasta Genes_By_Sample/"$PublicName"-pt_genes.fasta
		cp ../Data/"$Sample"_pt.fasta Contigs/"$PublicName"-pt_contigs.fasta
	done < ../$Label_file
	### Copy low and high-threshold trees and alignments
	mkdir -p Low_Thresold_Tree; mkdir -p High_Thresold_Tree;
	mkdir -p Low_Thresold_Tree/Aligned_Genes; mkdir -p High_Thresold_Tree/Aligned_Genes;
	mkdir -p Low_Thresold_Tree/Aligned_Genes_Trimmed; mkdir -p High_Thresold_Tree/Aligned_Genes_Trimmed;
	lt="$(echo $Label_file | sed 's/_DataDistribution.txt//g')"_0.1
	ht="$(echo $Label_file | sed 's/_DataDistribution.txt//g')"_clean_0.5
	echo $lt
	cp -a ../$lt/Aligned_genes/. Low_Thresold_Tree/Aligned_Genes/
	cp -a ../$lt/Aligned_genes/. Low_Thresold_Tree/Aligned_Genes_Trimmed/
	echo $ht
	cp -a ../$ht/Aligned_genes/. High_Thresold_Tree/Aligned_Genes/
	cp -a ../$ht/Aligned_genes/. High_Thresold_Tree/Aligned_Genes_Trimmed/
	
	rm AllSamples_Allgenes_Sel_SeqTable.csv
	rm AllSamples_Allgenes.fasta
	### Compress folder
	cd ..
	# tar -zcvf "$project"_dryad.tar.gz $project/
	zip -r "$project"_dryad.zip $project/
fi

echo "Done"





# http://www.iqtree.org/doc/Advanced-Tutorial
# iqtree -s concat.fasta -p concat.partitions -m GTR+F+I+G4 --prefix test_concat -T 4
#iqtree -p Aligned_genes_trimmed -m MFP --prefix test_p -B 1000 -T AUTO
# iqtree -s Concat_aligned_38genes.fasta -m MFP -B 1000 -T AUTO