#!/bin/bash
#SBATCH --job-name="AMAS"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=16000
#SBATCH --partition=medium
ncpu=4

#JHI
source activate phylo

dir=$1

AMAS.py summary -f fasta -d dna -i $dir/*.fasta -c 4
