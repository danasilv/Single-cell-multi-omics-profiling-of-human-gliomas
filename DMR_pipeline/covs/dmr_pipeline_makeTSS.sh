#!/usr/bin/env bash
#$ -cwd
#$ -N JOB_NAME
#$ -pe smp 1 -l mem=40G
#$ -l h_vmem=40G
#$ -m abe
#$ -M your_email
#$ -p -10

module load R/3.6.1
module load bedtools/2.27.1

### change path to scripts to folder where you will have the scripts
scripts=/scripts

### Intersect with UCSC +/- 1kb TSS
for i in *binarize.cov;
do
bedtools intersect -a ${i} -b $scripts/ucsc_refseq_V2.1kb.bed -wa -wb > ${i%.*}.promoters.1kb.cov;
done

### Create 2 Matrices

### A) A gene by cell methylation matrix 
Rscript $scripts/promoters_meth_by_cell.R ./ TSS_1kb

### B) A gene by cell number of CpG matrix 
Rscript $scripts/promoters_cpgs_by_cell.R ./ TSS_1kb