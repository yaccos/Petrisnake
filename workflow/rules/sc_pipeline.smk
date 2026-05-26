def get_lane_files_for_demultiplexing(wildcards):
    sample = wildcards.sample
    lanes = sample_lanes[sample]
    template = "{results}/{sample}/{sample}_QF_{lane}_R1.fastq"
    return [f"results/{sample}/{sample}_QF_{lane}_R1.fastq" for lane in lanes]


rule demultiplex:
    input:
        fastq=get_lane_files_for_demultiplexing,
        bc1=f"{barcode_dir}/BC1_5p_anchor_v2.fa",
        bc2=f"{barcode_dir}/BC2_anchored.fa",
        bc3=f"{barcode_dir}/BC3_anchored.fa",
    output:
        barcode_table="results/{sample}/{sample}_barcode_table.txt",
        bc_frame=temp("results/{sample}/{sample}_bc_frame.rds"),
        freq_table="results/{sample}/{sample}_frequency_table.txt",
    params:
        chunk_size=lambda wildcards: processed_config[wildcards.sample]["chunk_size"],
        adapter_sequence="GGTCCTTGGCTTCGC",
    log:
        "logs/{sample}/demultiplex.log",
    container:
        "docker://yaccos/posdemux:0.99.8"
    script:
        "../scripts/demultiplexer.R"


rule create_bc_plots:
    input:
        "results/{sample}/{sample}_frequency_table.txt",
    output:
        histogram="results/{sample}/{sample}_ReadsPerBC.pdf",
        knee_plot="results/{sample}/{sample}_kneePlot.pdf",
    params:
        bc_cutoff=lambda wildcards: processed_config[wildcards.sample]["bc_cutoff"],
    log:
        "logs/{sample}/bc_plots.log",
    container:
        "docker://yaccos/posdemux:0.99.8"
    script:
        "../scripts/create_bc_plots.R"


rule select_reads:
    input:
        freq_table="results/{sample}/{sample}_frequency_table.txt",
        barcode_table="results/{sample}/{sample}_barcode_table.txt",
    output:
        selected_barcode_table="results/{sample}/{sample}_selected_frequency_table.txt",
        barcode_database=temp("results/{sample}/{sample}_selected_barcode_table.sqlite"),
    log:
        "logs/{sample}/select_reads.log",
    params:
        bc_cutoff=lambda wildcards: processed_config[wildcards.sample]["bc_cutoff"],
        chunk_size=lambda wildcards: processed_config[wildcards.sample]["chunk_size"],
    container:
        "docker://yaccos/posdemux:0.99.8"
    script:
        "../scripts/select_reads_to_keep.R"
