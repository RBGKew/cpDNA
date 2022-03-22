# Plastid Dataset and Tree building Pipeline

Here we describe the pipeline used to produce the plastid dataset and tree published in Leempoel *et al.* (preprint).

The validated dataset is available on Dryad (link) and the annotated plastid genomes and contigs are available on INSDC under our umbrella project PRJEB35285

## Pipeline

### Input

A list of samples must be provided in a space delimited file, which should have the following format `[Project_Name]_labels.txt` and have two columns (without header): Sample_Name Label. 

Plastid contigs are expected to be located in the Data/ folder and Files should have the following format `[Sample_Name]_pt.fasta`

### Parameters

- *Reference*=References/cpDNA_RefCDS.fasta
- *min_cov*=50  Minimum coverage of the matching reference gene to keep a blast match
- *min_len_r*=0.2 minimum proportion of sum of genes length to keep a sample
- *min_sample_r*=0.5 minimum proportion of samples to keep a gene
- *slurmThrottle*=80 Number of parallel slurm jobs

### Running modes

The script has 3 running modes. In all cases, a label file must be provided

1. **filter**, to launch the gene recovery, followed by a filtering based on the minimum proportion of sum of genes. 

- Samples' contigs are blasted against the *Reference*. A custom python script `scripts/targets_from_blast.py` processes the blastn output by selecting the highest bitscore by gene and discarding matches below the minimum reference coverage  *min_cov*. An attempt is made a recovering all exons for chloroplast genes with multiple exons.
- All recovered gene sequences are concatenated into a single fasta file. The latter is then processed with `scripts/Transpose_s2g.py`, which writes fasta files per gene with all accepted samples. Here thresholds *min_len_r* & *min_sample_r* define the proportion of samples and genes kept.

2. **tree**, to build the plastid tree from recovered genes

- Genes are aligned with `mafft --adjustdirection --auto` and trimmed with `trimal  -gappyout`
- The tree is built with `iqtree -s Aligned_genes_trimmed/ -m GTR+G -B 1000` and labelled.

3. **distrib**, to assemble the validated dataset into a single compressed archive

Examples

```shell
cpDNA_pipeline.sh Araceae_labels.txt "tree" 0.1 #where Araceae_labels is the label file and 0.1 is the min_len_r, corresponding to 10% of the sum of genes length
cpDNA_pipeline.sh cpTree_v7_tree_DataDistribution.txt "distrib" 0 
```



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
