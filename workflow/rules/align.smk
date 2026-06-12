bwa_index_extensions = [".bwt", ".amb", ".pac", ".sa", ".ann"]


rule bwa_index:
    input:
        "{genome_dir}/{genome_basename}.{genome_ext}",
    output:
        temp(
            multiext(
                "{genome_dir}/{genome_basename}.{genome_ext}", *bwa_index_extensions
            )
        ),
    log:
        "logs/bwa_index/{genome_dir}/{genome_basename}.{genome_ext}_index.log",
    conda:
        "../envs/bwa.yml"
    shell:
        "bwa index {input} 2> {log}"


get_bwa_index = lambda wildcards: multiext(
    processed_config[wildcards.sample]["genome"], *bwa_index_extensions
)


# BWA alignment
rule bwa_align:
    input:
        reads="results/{sample}/{sample}_QF_{lane}_R2.fastq",
        genome=lambda wildcards: processed_config[wildcards.sample]["genome"],
        index=get_bwa_index,
    output:
        temp("results/{sample}/{sample}_bwa_{lane}.sam"),
    log:
        "logs/{sample}/bwa_mem_{lane}.log",
    # Ensures there is a core left to do demultiplexing
    threads: max(1, workflow.cores - 1)
    conda:
        "../envs/bwa.yml"
    shell:
        "bwa mem -a -k 12 -Y -t {threads} {input.genome}  {input.reads} > {output} 2> {log}"
