configfile: "config.yaml"

rule sites_only:
	input:
		"site_matrix.csv",
        "Pairwise_stats_rrbs.csv"

rule by_region:
	input:
		"region_matrix_pct.csv",
		"region_matrix_count.csv",
		"Pairwise_stats_rrbs.csv"

rule demultiplex:
	input:
		r1="{prefix}_R1_001.fastq",
		r2="{prefix}_R2_001.fastq"
	output:
		expand("{{prefix}}.{barcode}.R1.fastq", barcode=config["barcodes"]),
		expand("{{prefix}}.{barcode}.R2.fastq", barcode=config["barcodes"])
	shell:
		"bash {config[scripts]}/demultiplex.sh {config[scripts]}/splitFastqPair.pl {input.r1}"

rule bismark_align:
	input:
		r1="{prefix}.{barcode}.R1.fastq",
		r2="{prefix}.{barcode}.R2.fastq"
	output:
		"{prefix}.{barcode}.R1.fastq_bismark_bt2_pe.bam"
	shell:
		"bismark --multicore 4 -X 1000 --path_to_bowtie /nfs/sw/bowtie2/bowtie2-2.2.8/ --un --ambiguous {config[genome]} -1 {input.r1} -2 {input.r2}"

rule bismark_extract:
	input:
		"{prefix}.{barcode}.R1.fastq_bismark_bt2_pe.bam"
	output:
		"{prefix}.{barcode}.R1.fastq_bismark_bt2_pe.bismark.cov.gz",
		"CpG_context_{prefix}.{barcode}.R1.fastq_bismark_bt2_pe.txt"
	shell:
		"bismark_methylation_extractor --bedgraph --comprehensive {input}"

rule pdr:
	input:
		expand("CpG_context_{prefix}.{barcode}.R1.fastq_bismark_bt2_pe.txt", prefix=config["prefixes"], barcode=config["barcodes"])
	output:
		"PDR_rrbs.txt"
	shell:
		"python {config[scripts]}/compute_pdr.py {input} > {output}"

rule qc_stats:
	input:
		"PDR_rrbs.txt"
	params:
		cells=expand("{prefix}.{barcode}", prefix=config["prefixes"], barcode=config["barcodes"])
	output:
		"QC_stats_rrbs.csv"
	shell:
		"python {config[scripts]}/cov_batch.py {params.cells} {input} rrbs"

rule pairwise_stats:
	input:
		"QC_stats_rrbs.csv"
	output:
		"Pairwise_stats_rrbs.csv"
	shell:
		"python {config[scripts]}/compute_pair.py {input} > {output}"

rule site_matrix:
	input:
		expand("{prefix}.{barcode}.R1.fastq_bismark_bt2_pe.bismark.cov.gz", prefix=config["prefixes"], barcode=config["barcodes"])
	output:
		"site_matrix.csv"
	shell:
		"python {config[scripts]}/build_features.py {input} > {output}"

rule site_bed:
	input:
		"site_matrix.csv"
	output:
		"site_matrix.bed"
	shell:
		"python {config[scripts]}/csv_to_site_bed.py {input} > {output}"

rule site_intersect:
	input:
		"site_matrix.bed"
	output:
		"site_intersect.bed"
	shell:
		"intersectBed -a {config[regions]} -b {input} -wa -wb > {output}"

rule region_matrix:
	input:
		"site_matrix.csv",
		"site_intersect.bed"
	output:
		"region_matrix_pct.csv",
		"region_matrix_count.csv"
	shell:
		"python {config[scripts]}/met_matrix.py {config[path]}"
