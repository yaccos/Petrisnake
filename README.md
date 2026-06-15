# Petrisnake: A secondary analysis pipeline for PETRI-seq data.

This is a Snakemake pipeline for the secondary computational analysis of single cell RNA-seq data from the PETRI-seq protocol (https://www.nature.com/articles/s41564-020-0729-6 and https://www.nature.com/articles/s41586-024-08124-2), this is: From the input FASTQ files, this workflow constructs a gene count table showing the expression of each gene in each cell. Petrisnake is available on WorkflowHub (https://workflowhub.eu/workflows/2081).

![Workflow rulegraph](./images/rulegraph.svg)

# Dependencies

This pipeline is only tested under Linux running on x86-64 architecture, although it might work on macOS and other CPU architectures if some work is put into dependency management. The workflow will complain, but still work if the system dependency `fuse2fs` is not installed. For running the pipeline the following dependencies are required (tested versions):

* [Snakemake](https://snakemake.github.io/) (9.14.5)
* [conda](https://www.anaconda.com/docs/getting-started/miniconda/main) (25.7.0)
* [Apptainer](https://apptainer.org/) (1.4.5)

The configuration tool GUI is a R Shiny application with the following dependencies:

* R (4.5.2)
* shiny (1.12.1)
* shinyBS (0.61.1)
* purrr (1.2.0)
* yaml (2.3.12)

For convenience, the file `config/environment.yml` describes the conda dependencies for the configuration GUI. In order to install it, run `conda env create --file config/environment.yml`. This should create a conda environment named `petrisnake-config` with the required dependencies.


The `posDemux` R package (https://github.com/yaccos/posDemux) is a keystone dependency of this workflow. This package is scheduled to be released as part of Bioconductor 3.23. For the pipeline itself, the package is available as a container image on Docker Hub (https://hub.docker.com/repository/docker/yaccos/posdemux/general) and is automatically used by Snakemake in the workflow provided that the argument `--software-deployment-method apptainer` is used.


When running the container, you may get the warning messages:
```
During startup - Warning messages:
1: Setting LC_CTYPE failed, using "C" 
2: Setting LC_COLLATE failed, using "C" 
3: Setting LC_TIME failed, using "C" 
4: Setting LC_MESSAGES failed, using "C" 
5: Setting LC_MONETARY failed, using "C" 
6: Setting LC_PAPER failed, using "C" 
7: Setting LC_MEASUREMENT failed, using "C"
```

These warnings can be safely ignored, but if you still want to silence them, run `export LC_ALL=C` in your shell before running Snakemake.

For using the cutoff selection tool, the `posDemux` package must be run outside the workflow. At the time of writing Bioconductor release 3.23 is not out yet and the package requires a minimum R version of 4.6.0 as per Bioconductor's policies. This makes the package trickier to install, but there are still some options:

* Install the devel version of R which is 4.6. You may use the `rig` version manager (https://github.com/r-lib/rig) for this as it allows you to install multiple versions of R and switch between them as needed.
* Download the sources from GitHub and manually edit the minimum version requirements to your installed R version (everything above 4.4 should work).
* Pull the container image with `posDemux` with Docker and use the bundled RStudio Server instance: `docker run -e PASSWORD=<YOUR_PASS> -p 8787:8787 yaccos/posdemux:0.99.8`. Then open `http://localhost:8787/` in your web browser, type `rstudio` as the username and `<YOUR_PASS>` as the password.

For the rules inside the workflow which do not use R scripts, Snakemakes automatically handles the conda environments for these rules such that no manual intervention is required for the dependency management. These conda environment files are located in `workflow/envs`.

# Configuration

The pipeline reads `config/config.yaml` and normalizes the values through `workflow/rules/handle_config.smk`. Every run must define the keys listed below; optional fields default to empty strings or can be overridden per sample when needed. The bundled Shiny helper (`config/config_shiny.R`) generates exactly this structure.

**Global keys**
- `prefix`: Base directory prepended to any relative path (for example `resources/`). Leave empty when using absolute paths. This option is special as it is not overridden by the key `prefix` inside the sample. Instead, the global and sample-specific prefixes are concatenated. In addition, this prefix relates to **all** file paths, not just the FASTQ files, but also the reference genome and reference annotation. 
- `reference_genome`: FASTA file relative to `prefix` unless the value is absolute. The file must exist in a writable directory so `bwa` can index it.
- `reference_annotation`: Annotation GTF/GFF file relative to `prefix` unless absolute.
- `suffix`, `forward_suffix`, `reverse_suffix`: Optional strings appended to every FASTQ path. They allow you to express shared filename parts one time (for example `_001.fastq.gz`).
- `streaming_chunk_size`: Integer, number of reads to load in each chunk while streaming FASTQ pairs during demultiplexing.
- `bc_cutoff`: Optional global barcode cutoff. When omitted you can still run `determine_bc_cutoff`; the value becomes mandatory before running `all`.
- `feature_tag`, `gene_id_attribute`: Required feature and identifier attribute names extracted from the reference annotation (for example `Coding_or_RNA` and `name`).
- `samples`: Mapping of sample identifiers to per-sample configuration dictionaries (see below). The mapping itself is mandatory.

**Sample entries**
- `prefix`: Prepended to the FASTQ file names in addition to the global `prefix`. Useful for subdirectories such as `random20000/random20000_S1_`.
- `suffix`, `forward_suffix`, `reverse_suffix`: Optional overrides. When omitted, the global suffixes are used.
- `reference_genome`, `reference_annotation`, `feature_tag`, `gene_id_attribute`, `bc_cutoff`, `streaming_chunk_size`: Optional overrides. Provide only what differs from the global defaults.
- `lanes`: Dictionary mapping lane labels to the strings found in the filenames (for example `L1: L001`). Lane keys must match the regex `^[^_]+$` (underscores are rejected). The workflow derives each FASTQ pair as `<prefix><forward_prefix><lane_id><forward_suffix><suffix>` and `<prefix><reverse_prefix><lane_id><reverse_suffix><suffix>`.

The minimal viable YAML therefore looks like:

```yaml
prefix: resources/
reference_genome: U00096_JE2.fa
reference_annotation: U00096_JE2_rRNA.gff
feature_tag: Coding_or_RNA
gene_id_attribute: name
streaming_chunk_size: 4000
samples:
	random20000:
		prefix: random20000/random20000_S1_
		lanes:
			L1: L001
		forward_suffix: _R1
		reverse_suffix: _R2
		suffix: _001.fastq.gz
		bc_cutoff: 3000
```

Each sample entry is validated at runtime. Missing required fields raise descriptive `ValueError`s (for example when a genome is neither defined globally nor for the sample). The `lanes` dictionary may contain multiple entries when samples are split across lanes. In this case, the workflow processes the FASTQ files for the forward reads for a sample seperately until demultiplexing where the files are streamed sequentially within a single rule execution. Likewise, the reverse read FASTQ files of a sample are processed separately until the results are combined at the final step of the workflow.

To create or update the configuration interactively, activate the GUI environment (`conda activate petrisnake-config`) and launch `Rscript config/config_shiny.R`. The app mirrors all of the fields documented above and writes a ready-to-run YAML file.

# Running the pipeline

Due to the fact that Snakemake uses Apptainer and Conda to handle the dependencies of the rules, the flag `--software-deployment-method apptainer conda` must be provided when launching the pipeline.

When the barcode threshold to use is unknown, we first run the pipeline to the end of the demultiplexer. In that case, run
`snakemake --software-deployment-method apptainer conda --cores <NUMBER_OF_CORES> determine_bc_cutoff`.

This will perform the demultiplexing of the forward reads and produce the following files in interest:

* The barcode table at `results/<sample>/<sample>_barcode_table.txt` which contains the barcode and UMI assignments to each of the reads.
* The frequency table at `results/<sample>/<sample>_frequency_table.txt` which lists the frequency of each barcode combination in decending order (most common ones on top). This table is used by `posDemux::interactive_bc_cutoff()` to determine the barcode cutoff, this is how many of the top barcode combinations are kept.
* The log file of the demultiplexing at `results/<sample>/demultiplex.log`. This file serves both as a progress indicator of the process and provides an informative summary of the results in the end.

After the barcode cutoff is determined and entered into the config file, we are ready to run the rest of the pipeline:

`snakemake --software-deployment-method apptainer conda --cores <NUMBER_OF_CORES> all`

The additional files of interest produced after the entire pipeline has finished are:

* The selected frequency table at `results/<sample>/<sample>_selected_frequency_table.txt`. This is a truncated and renormalized version of the frequency table containing only the reads within the barcode cutoff.
* The figures of the knee plot and density of barcode frequencies located at `results/<sample>/<sample>_kneePlot.pdf` and `results/<sample>/<sample>_ReadsPerBC.pdf`. The dashed vertical line indicates the selected cutoff. Note that the density plot is scaled by the number of reads on the y-axis in order to make it easier to interpret. In case you are unhappy with how the plots are made, you can recreate the plots from the frequency table using `posDemux::knee_plot()` and `posDemux::freq_plot()`.
* The gene count matrix at `results/<sample>/<sample>_gene_count_matrix.txt`. This is a tab-seperated file (remember to explicitly set `\t` as the delimiter) showing the number of detected UMIs per gene for all cells. The genes form the columns, whereas the cells form the rows.
* The UMI count table at `results/<sample>/<sample>_umi_count_table.txt` which can be seen as a linear form of the gene count matrix. Each row contains the cell barcode, the UMI, its associated gene and how many duplicates were found of this UMI. 
* The log file from the gene counting at `results/<sample>/count_genes.log` which serves both as a progress indicator and provides a summary of the gene count table.

This repository is already prepared with an example dataset `random20000` and a corresponding pre-populated config file, meaning that the user can try the pipeline right out of the box.
