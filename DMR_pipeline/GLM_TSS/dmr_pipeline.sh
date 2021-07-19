#!/usr/bin/env bash
#$ -cwd
#$ -N JOB_NAME
#$ -pe smp 1 -l mem=5G
#$ -l h_vmem=5G
#$ -m abe
#$ -M your_email
#$ -p -10

module load R/3.6.1
module load bedtools/2.27.1

### change path to scripts to folder where you will have the scripts
scripts=/scripts

### Perform DMR analysis
Rscript allcells_promoter_per_samp_GLM_global.R TSS_1kb_binarize_promoters_meth.rds TSS_1kb_binarize_promoters_cpgs.rds Diff Stem