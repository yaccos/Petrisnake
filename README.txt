Benchmark scripts for the paper "An improved computational pipeline for single-cell RNA-seq of bacteria",
 benchmarking the old PETRI-seq pipeline by Blattman et al.

Original pipeline written by Sydney Blattman
Modifications carried out by Jakob Peder Pettersen:
* Adapted the pipeline to use Python 3 instead of Python 2.7
* Made cosmetical changes to the pipeline scripts to make it more readable, utilizing f-strings in Python. The functionality should be the same.
* Changed the dataset being used to SRR28148450 which stems from the paper by Yan et al. 2024
* Changed featureCounts_directional_5.py to use the same reference annotation settings as Yan et al.
* Changed the aligner from bwa aln to bwa mem
Last updated: 02-06-2026

Dependencies:
python 3.12.12
fastqc (v0.12.1)
cutadapt (v5.2)
seqkit v2.11.0
GNU parallel (20251122)
UMI-tools (v1.1.6)
subread (v2.1.1, contains featureCounts)
bwa (0.7.19)
samtools (1.7)
PEAR (0.9.6)
pandas (2.3.3)
matplotlib (3.10.8)
numpy (2.3.5)
GNU time (1.9)

How to run the benchmark:

1. Gather the dependencies:
 The quickest and easiest way to reproduce the environment is to use conda.
 Run: conda env create --file environment.yml
 This should create a conda environment named blattman-petriseq containing the dependencies for the pipeline
2. Download the sequence data: Due to the sheer size of the dataset,
 the sequence files are not included in this repo. Download the sequence data for SRR28148450 from NCBI and put the resulting files
 SRR28148450_1.fastq and SRR28148450_2.fastq in benchmark/EC_3h_deep
3. Navigate to the benchmark directory
4. If using conda, run: conda activate blattman-petriseq
5. Run the benchmark, starting the benchmark by running: command time -v -o benchmark_timing.log ./run_benchmark.sh
6. The benchmark timings should be available as benchmark_timing.log, the gene count table as EC_3h_deep_v11_threshold_0_mixed_species_gene_matrix.txt, and the
 UMI summary table as EC_3h_deep_v11_threshold_0_filtered_mapped_UMIs.txt
