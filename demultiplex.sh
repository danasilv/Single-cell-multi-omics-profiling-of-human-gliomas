#!/bin/bash

DEMULT=$1
INPUT=$2
BASENAME=`basename $INPUT _R1_001.fastq`
echo $BASENAME

perl $DEMULT ${BASENAME}_R1_001.fastq ${BASENAME}_R2_001.fastq

for file in *${BASENAME}*.fastq.*.fastq; do
    mv "$file" $(echo "$file" | sed -E 's/(.*)_(R1|R2)_001\.fastq\.(.*)\.fastq/\1\.\3\.\2\.fastq/')
done
