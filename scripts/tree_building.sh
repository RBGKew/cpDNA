#!/bin/bash

#SBATCH --job-name="make_tree"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=20
#SBATCH --ntasks=1
#SBATCH --mem=120000		### Paul B. - increased from 80GB to 120 GB

## Kew
# module load python/3.7.9
# module load amas/1.0
# module load iqtree/2.0.6
# module load fasttree
# module load trimal/1.4.1
# module load mafft/7.471
## JHI
### source activate phylo	# Paul B. - changed to set env outside of this script

project=$1
Label_file=$2
ncpu=28			# Paul B. - increased from 20 to 28 cpu


# fasttree -gtr -nt -gamma -log "$project"_fasttree.log < concat.fasta > "$project"_fasttree.nwk
# nw_rename -l "$project"_fasttree.nwk ../$Label_file > "$project"_fasttree_labelled.nwk
date +"%T"
# Paul B. changed from iqtree to use iqtree2, the name of the native program:
iqtree2 -s Aligned_genes_trimmed/ -m GTR+G -B 1000 -T $ncpu --prefix "$project" --seqtype DNA -v
date +"%T"
nw_rename -l "$project".treefile ../$Label_file > "$project"_iqtree_labelled.nwk