#!/usr/bin/env bash

META_PATH=/path_to_RNA_FASTQs/

# Note: each directory specified should have a sub-directory called 'fastq' containing the fastq files for the sample

for directory in ${META_PATH}/*; do
    echo $directory
    qsub run_RNA.sh $directory
done
