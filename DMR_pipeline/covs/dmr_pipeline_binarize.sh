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

### 1) Perform 50k CpG filtration (this can be skipped if you already place in the covs folder the cells with >50k CpGs)
### covs with <50k are placed in "failed_QC"

bash filter_cpgs.sh

### 2) Binarize cov files

Rscript $scripts/binarize_cov.R ./