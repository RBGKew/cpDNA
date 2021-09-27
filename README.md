# Chloroplast Tree Of Life Pipeline

Here we describe the pipeline used to produce the chloroplast tree from Leempoel *et al.* (preprint). The pipeline was used both for the pre and final tree.

## Input

Chloroplast contigs are located in the Data/ folder. Files should have the following format `Sample_pt.fasta`

Samples for which chloroplast genes are already extracted (skipping blastn to identify chloroplast genes) are located in cpGenes_by_sample/. Files should have the following format `Sample_cpGenes.fasta`

## Parameters

- *Reference*=References/cpDNA_RefCDS.fasta
- *min_cov*=50  Minimum coverage of the matching reference gene to keep a blast match
- *min_len_r*=0.2 minimum proportion of sum of genes length to keep a sample
- *min_sample_r*=0.5 minimum proportion of samples to keep a gene
- *slurmThrottle*=80 Number of parallel slurm jobs

## Dependencies

```python
- blastn
- iqtree/2.1.3
- mafft/7.487
- trimal/1.3
- newick-utils 1.6 
- python3.6
  - Pandas
  - Numpy
  - Bio
```

## Pipeline

Pipeline is launched with a label file, which should have the following format `project_labels.txt` and has two columns (without header) separated by a space: sample label

```shell
sbatch cpDNA_pipeline.sh cpTree_v3_pretree_labels.txt
```

The pipeline includes the following steps:

1. Samples' contigs are blasted against the *Reference*. A custom python script `scripts/targets_from_blast.py` processes the blastn output by selecting the highest bitscore by gene and discarding matches below the minimum reference coverage  *min_cov*. An attempt is made a recovering all exons for chloroplast genes with multiple exons.
2. All recovered gene sequences are concatenated into a single fasta file. The latter is then processed with `scripts/Transpose_s2g.py`, which writes fasta files per gene with all accepted samples. Here thresholds *min_len_r* & *min_sample_r* define the proportion of samples and genes kept.
3. Genes are aligned with `mafft --adjustdirection --auto` and trimmed with `trimal  -gappyout`
4. The tree is built with `iqtree -s Aligned_genes_trimmed/ -m GTR+G -B 1000` and labelled.
