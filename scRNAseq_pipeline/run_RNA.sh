#!/usr/bin/env bash
# -V
# -cwd
# -l h_rss=5G

module load star/2.5.2a
module load R/3.6.1

MY_PATH=$(readlink -f ${1})
echo $MY_PATH
name="${MY_PATH##*/}"
echo $name
cd ${MY_PATH}/fastq

echo "path: ${MY_PATH}/fastq" > config.yaml
echo "name: ${name}" >> config.yaml
echo "samples:" >> config.yaml

for file in *R1.fastq; do
    sample=`basename $file _R1.fastq`
    echo "- ${sample}" >> config.yaml
done

MY_SNAKEMAKE=$(which snakemake)

$MY_SNAKEMAKE -s RNA_pipeline_STAR.mk --configfile ${MY_PATH}/fastq/config.yaml --cluster-config cluster.yaml --cluster "qsub -V -cwd -N hg38_RNA -pe smp 1 -l h_rt={cluster.rt},h_rss={cluster.mem} -p -10 -o qsub_output -e qsub_error" --jobs 1000
