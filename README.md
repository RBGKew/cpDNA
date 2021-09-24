# Chloroplast Tree Of Life Pipeline

Here we describe the pipeline used to produce the chloroplast tree from Leempoel et al 2022.

## Input

Chloroplast contigs are located in the Data/ folder. Files should have the following format `Sample_pt.fasta`

Samples for which chloroplast genes are already extracted (skipping blastn to identify chloroplast genes) are located in cpGenes_by_sample/. Files should have the following format `Sample_cpGenes.fasta`

## Parameters

Reference=References/cpDNA_RefCDS.fasta
min_cov=50  Minimum coverage of a reference gene to keep a blast match
min_gene_r=0.8 Quantile of sum of gene length (e.g. 0.8 means the minimum sum of gene length will be 20% of the top 20% percentile)
min_len_r=0.2 minimum % of sum of gene length calculated with min_gene_r
min_sample_r=0.5 minimum % of samples to keep a gene
slurmThrottle=80 Number of parallel jobs run

## Dependencies

```python
- blastn
- iqtree 2
- trimAI
- python3
  - Pandas
  - Numpy
  - Bio
```

## Pipeline

Pipeline is launched with a label file, which should have the following format `project_labels.txt`

```shell
sbatch cpDNA_pipeline.sh cpTree_v3_pretree_labels.txt
```

