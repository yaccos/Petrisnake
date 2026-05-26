log_file = snakemake@log[[1]]
log_handle  <- file(log_file, open = "w")
sink(log_handle, append = TRUE, type = "output")
sink(log_handle, append = TRUE, type = "message")

suppressMessages({
    library(posDemux)
    library(Biostrings)
    library(purrr)
    library(glue)
    library(tibble)
    library(ShortRead)
})

BARCODE_WIDTH <- 7L
ALLOWED_MISMATCHES <- 1L
DEFAULT_BC3_TO_BC2_ADAPTER <- "GGTCCTTGGCTTCGC"

log_progress <- function(msg) {
    message(glue("{date()} => {msg}"))
}

filter_fastq_for_adapter <- function(filepath, adapter_sequence, chunk_size) {
    log_progress(glue("Filtering {filepath} for adapter {adapter_sequence}"))

    output_file <- tempfile(
        pattern = paste0(tools::file_path_sans_ext(basename(filepath)), "_adapter_"),
        fileext = ".fastq"
    )
    streamer <- ShortRead::FastqStreamer(filepath, n = as.integer(chunk_size))
    on.exit(close(streamer), add = TRUE)

    wrote_any_reads <- FALSE
    total_reads <- 0L
    retained_reads <- 0L

    repeat {
        chunk <- ShortRead::yield(streamer)
        if (length(chunk) == 0L) {
            break
        }

        total_reads <- total_reads + length(chunk)
        keep_reads <- Biostrings::vcountPattern(
            adapter_sequence,
            ShortRead::sread(chunk),
            fixed = TRUE
        ) > 0L
        kept_in_chunk <- sum(keep_reads)
        retained_reads <- retained_reads + kept_in_chunk

        if (kept_in_chunk > 0L) {
            ShortRead::writeFastq(
                chunk[keep_reads],
                file = output_file,
                mode = if (wrote_any_reads) "a" else "w"
            )
            wrote_any_reads <- TRUE
        }
    }

    if (!wrote_any_reads) {
        file.create(output_file)
    }

    log_progress(glue(
        "Retained {retained_reads} of {total_reads} reads from {filepath}"
    ))

    output_file
}

bc_frame <- tibble(bc_name = glue("bc{1:3}"))
bc_frame$filename <- snakemake@input[bc_frame$bc_name] |> as.vector()

input_file <- snakemake@input[["fastq"]] |> as.vector()

output_table_file <- snakemake@output[["barcode_table"]]
output_bc_frame  <- snakemake@output[["bc_frame"]]
output_freq_table  <- snakemake@output[["freq_table"]]
chunk_size <- snakemake@params[["chunk_size"]]
adapter_sequence <- snakemake@params[["adapter_sequence"]]
if (is.null(adapter_sequence) || !nzchar(adapter_sequence)) {
    adapter_sequence <- DEFAULT_BC3_TO_BC2_ADAPTER
}
adapter_sequence <- toupper(adapter_sequence)

sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")

segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)

bc3_to_bc2_adapter <- "GGTCCTTGGCTTCGC"

bc_frame$stringset <- map(bc_frame$filename, function(filepath) {
    glue("Reading filepath {filepath}")  |>  log_progress()
    raw_stringset <- Biostrings::readDNAStringSet(filepath = filepath)
    glue("Trimming away adapters in barcode file")  |> log_progress()
    # The FASTA files contain the barcodes in addition to the adapters, so we must filter them out
    Biostrings::subseq(raw_stringset, start = 1L, width = BARCODE_WIDTH)
})

names(bc_frame$stringset) <- bc_frame$bc_name

filtered_input_file <- map_chr(
    input_file,
    filter_fastq_for_adapter,
    adapter_sequence = adapter_sequence,
    chunk_size = chunk_size
)

callbacks <- streaming_callbacks(input_file = filtered_input_file,
                                 output_table_file = output_table_file,
                                 chunk_size = chunk_size,
                                 verbose = TRUE)

streaming_res <- rlang::exec(streaming_demultiplex, !!! callbacks,
                             barcodes=bc_frame$stringset |> rev(), allowed_mismatches = ALLOWED_MISMATCHES,
            segments = sequence_annotation, segment_lengths = segment_lengths)

freq_table  <- streaming_res$freq_table
log_progress("Writing frequency table...")
write.table(x = freq_table, file = output_freq_table, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
log_progress("DONE")
cat("\n")
print(streaming_res$summary_res)

saveRDS(object = bc_frame, file = output_bc_frame, compress = FALSE)

sink(type="message")
sink(type="output")
close(log_handle)
