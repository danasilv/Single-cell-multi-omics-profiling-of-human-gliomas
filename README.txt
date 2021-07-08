README.txt

Dependencies:
- Local installation of snakemake (see website for instructions)
- python
- bismark/0.14.5
- bedtools/2.25.0
- samtools/1.3.1

Execution:
bash ProcessRRBSLauncher.sh --argument value --flag

Required arguments:
-f/--fastq: Path to directory containing unzipped fastq files of the following format: <pool_prefix>_R1_001.fastq and <pool_prefix>_R2_001.fastq
-g/--genome: Path to bismark-compatible genome
-p/--pipeline: Path to directory containing pipeline (ProcessRRBS.mk) and other scripts
-t/--target: Target of pipeline (See potential values below)

Potential Values of target argument (in reverse order of execution):
- by_region: execute full pipeline and create region by cell matrix (requires optional regions argument below)
- sites_only: execute full pipeline up to site by cell matrix creation
- bismark_extract: execute pipeline up to bismark methylation extraction
- bismark_align: execute pipeline up to bismark alignment
- demultiplex: execute pipeline up to demultiplexing of fastq files

Optional arguments:
-r/--regions: Bed file containing genomic regions to analyze (required if target specified as 'by_region')
-c/--cluster-conifg: path to config file if submitting to cluster (see examples in this directory)

Optional flags:
--qsub: submit underlying jobs to cluster (requires cluster configuration file provided to -c argument)
--dry-run: map out pipeline but do not execute
--dag: create directed acyclic graph representing pipeline