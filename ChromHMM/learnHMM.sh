#!/usr/bin/bash
#$ -cwd
#$ -N JOB_NAME
#$ -pe smp 1 -l mem=20G
#$ -l h_vmem=20G
#$ -m abe
#$ -M your_email

java -mx1200M -jar Spectacle.jar LearnModel -chromhmm /BINARIZE/ /HMM 18 hg38