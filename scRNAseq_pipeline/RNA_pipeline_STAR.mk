rule all:
    input:
       expand("{name}.tpm.counts", name=config["name"])

rule RSEM:
    input:
        R1="{sample}_R1_001.fastq",
        R2="{sample}_R2_001.fastq"
    output:
        "{sample}.RSEM.genes.results"
    shell:
        "/RSEM/RSEM-1.3.1/rsem-calculate-expression \
                        --paired-end --star --star-path /path_to_star \
                        -p 1 {config[path]}/{input.R1} {config[path]}/{input.R2} \
                        /path_to_ref_genome \
                        {config[path]}/{wildcards.sample}.RSEM"

rule generate_matrix:
    input:
        file=expand("{dataset}.RSEM.genes.results", dataset=config["samples"])
    output:
        matrix=expand("{name}.tpm.counts", name=config["name"])
    script:
        "build_matrix_STAR_tpm.R"