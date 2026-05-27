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

log_progress <- function(msg) {
    message(glue("{date()} => {msg}"))
}

bc_frame <- tibble(bc_name = glue("bc{1:3}"))
bc_frame$filename <- snakemake@input[bc_frame$bc_name] |> as.vector()

input_file <- snakemake@input[["fastq"]]

output_table_file <- snakemake@output[["barcode_table"]]
output_bc_frame  <- snakemake@output[["bc_frame"]]
output_freq_table  <- snakemake@output[["freq_table"]]
chunk_size <- snakemake@params[["chunk_size"]]

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

loader <- function(state) {
    if (!state$output_table_initialized) {
        message("Starting to load sequences")
    }
    chunk <- Biostrings::readDNAStringSet(
        filepath = input_file,
        format = "fastq",
        nrec = chunk_size,
        skip = state$total_reads
    )
    n_reads_in_chunk <- length(chunk)


    if (n_reads_in_chunk == 0L && state$output_table_initialized) {
        # The case when the initial chunk is empty is given
        # special treatment since
        # we want to create the table regardless
        # No more reads to process, make the outer framework return
        final_res <- list(
            state = state,
            sequences = NULL,
            should_terminate = TRUE
        )
        message("Done demultiplexing")
        return(final_res)
    }
    state$total_reads <- state$total_reads + n_reads_in_chunk
    chunk_without_adapter <- Biostrings::vcountPattern(bc3_to_bc2_adapter, chunk, max.mismatch = 1L) |> (\(x) chunk[x == 0L])()
    state$reads_without_adapter  <- state$reads_without_adapter + length(chunk_without_adapter)
    list(
        state = state,
        sequences = chunk,
        should_terminate = FALSE
    )
}

state_init <- list(
    total_reads = 0L, reads_without_adapter = 0L, demultiplexed_reads = 0L,
    output_table_initialized = FALSE
)

archiver <- function(state, filtered_res) {
    barcode_matrix <- filtered_res$demultiplex_res$assigned_barcodes
    barcode_names <- colnames(barcode_matrix)
    read_names <- rownames(barcode_matrix)
    # If the table has no rows, we may risk getting a NULL value
    if (is.null(read_names)) {
        read_names <- character()
    }
    barcode_table <- as.data.frame(barcode_matrix)
    read_name_table <- data.frame(read = read_names)

    chunk_table <- cbind(read_name_table, barcode_table)
    if (!state$output_table_initialized) {
        append <- FALSE
        state$output_table_initialized <- TRUE
    } else {
        append <- TRUE
    }

    readr::write_tsv(
        x = chunk_table,
        file = custom_output_table,
        append = append,
        col_names = !append,
        eol = "\n"
    )
    state <- within(state, {
        demultiplexed_reads <- demultiplexed_reads + nrow(barcode_matrix)
        paste0(
            "Processed {total_reads} reads,", " ",
            "successfully demultiplexed {demultiplexed_reads} reads so far..."
        ) %>%
            glue::glue() %>%
            message()
    })
    state
}



callbacks <- streaming_callbacks(input_file = input_file,
                                 output_table_file = output_table_file,
                                 chunk_size = chunk_size,
                                 verbose = TRUE)

streaming_res <- rlang::exec(streaming_demultiplex, state_init = state_init,
                            loader = loader, archiver = archiver,
                             barcodes=bc_frame$stringset |> rev(), allowed_mismatches = 1L,
            segments = sequence_annotation, segment_lengths = segment_lengths)

final_state  <- streaming_res$state_final
total_reads  <- final_state$total_reads
reads_without_adapter <- final_state$reads_without_adapter
lacking adapter_percentage  <- round(total_reads / reads_without_adapter * 100, 2L)
demultiplexed_reads  <- final_state$demultiplexed_reads
demultiplexing_success  <- round(demultiplexed_reads / reads_without_adapter * 100, 2L)


glue(
    "Read of total of {total_reads} reads,
    of which {reads_without_adapter} ({lacking_adapter_percentage}%) lacked adapter"
    )

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
