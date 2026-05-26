# Kinetics-seq
Kinetics-seq enables comprehensive profiling of single-cell RNA kinetics in vivo to reveal dynamic tumor heterogeneity

## Abstract
Kinetics-seq is a temporal analysis technique combining in vivo metabolic labeling with single-cell RNA sequencing. It quantifies dynamic parameters including RNA synthesis and degradation at single-cell resolution, overcoming the limitation of conventional scRNA-seq that only captures static gene expression. This method helps dissect transcriptional heterogeneity and gene regulatory patterns in tumor cells, and screens active genes and pathways closely linked to tumor progression. It serves as a powerful tool for tumor mechanism research, biomarker identification and targeted therapy development.

## Set Up Environment
```
System: centOS
python 3.8.18
R 4.2.1
samtools 1.6
ngs-tools 1.8.5
pysam 0.20.0
```

## Obtaining the Expression Matrix for Kinetics-seq
1. 01_Preprocessing: Raw FASTQ files were processed using DynamicEX to generate the initial expression matrix. Aligned BAM files were mapped to the standard reference genome to quantify the expression matrices of newly transcribed and pre-existing RNAs. A binomial mixture model was applied to correct the expression signals of new RNAs from Kinetics-seq data and estimate the substitution rate induced by metabolic RNA labeling.
2. 02_computing: Kinetic parameters including RNA synthesis and degradation rates of each gene in individual cells were calculated based on the newly transcribed RNA matrix.
3. 03_analysis: Scripts for generating all figures presented in the manuscript.


