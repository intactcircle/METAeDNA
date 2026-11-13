#' Run DADA2 on paired FASTQ files
#'
#' Filter paired reads, learn error models, denoise both directions, merge the
#' pairs, and remove chimeras. The input may come from `barcode_split()` or
#' from any other source that provides paired FASTQ files by sample.
#'
#' @param forward_reads Character vector of forward FASTQ paths.
#' @param reverse_reads Character vector of reverse FASTQ paths in the same
#'   sample order as `forward_reads`.
#' @param output_dir Empty directory in which DADA2 results are written.
#' @param sample_names Sample names. By default, names are derived from the
#'   forward FASTQ filenames.
#' @param trunc_len R1 and R2 truncation lengths. `c(0, 0)` keeps full reads.
#' @param max_ee Maximum expected errors for R1 and R2.
#' @param n_threads Number of threads, or `TRUE` to let DADA2 choose.
#' @return A list containing the final `sequence_table` and read `tracking`.
#'   It is also saved to `dada2_result.rds`.
#'
#' @examples
#' \dontrun{
#' barcode_split(
#'   "run_R1.fastq.gz", "run_R2.fastq.gz", "samples.csv", "split_fastq"
#' )
#' forward_reads <- list.files("split_fastq", "_R1[.]fastq[.]gz$",
#'   full.names = TRUE)
#' reverse_reads <- list.files("split_fastq", "_R2[.]fastq[.]gz$",
#'   full.names = TRUE)
#'
#' dada2_result <- dada2_pipeline(
#'   forward_reads,
#'   reverse_reads,
#'   output_dir = "dada2_results",
#'   n_threads = 8
#' )
#' dada2_result$sequence_table
#' }
#' @export
dada2_pipeline <- function(forward_reads, reverse_reads, output_dir,
                           sample_names = NULL, trunc_len = c(0L, 0L),
                           max_ee = c(2, 2), n_threads = TRUE) {
  started <- Sys.time()
  on.exit(show_time("", 1, as.numeric(difftime(Sys.time(), started, units = "secs"))), add = TRUE)
  if (is.logical(n_threads)) .check_flag(n_threads, "n_threads") else .check_number(n_threads, "n_threads", 1, .Machine$integer.max, TRUE)
  .check_path(output_dir, "output_dir")
  if (!is.character(forward_reads) || !is.character(reverse_reads) || anyNA(c(forward_reads, reverse_reads))) {
    stop("Read paths must be character vectors without missing values.", call. = FALSE)
  }
  if (!requireNamespace("dada2", quietly = TRUE)) {
    stop("Install DADA2 with BiocManager::install('dada2').", call. = FALSE)
  }
  if (length(forward_reads) != length(reverse_reads) || !length(forward_reads)) {
    stop("forward_reads and reverse_reads must have the same non-zero length.",
         call. = FALSE)
  }
  input_files <- c(forward_reads, reverse_reads)
  if (any(!file.exists(input_files))) {
    stop("Input file not found: ",
         paste(input_files[!file.exists(input_files)], collapse = ", "),
         call. = FALSE)
  }
  if (length(trunc_len) != 2L || length(max_ee) != 2L) {
    stop("trunc_len and max_ee must each contain an R1 and R2 value.",
         call. = FALSE)
  }
  for (value in trunc_len) .check_number(value, "trunc_len", 0, .Machine$integer.max, TRUE)
  for (value in max_ee) .check_number(value, "max_ee", 0)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                  no.. = TRUE))) {
    stop("output_dir must be empty.", call. = FALSE)
  }

  if (is.null(sample_names)) {
    sample_names <- sub("_R1(?:\\.fastq|\\.fq)(?:\\.gz)?$", "",
                        basename(forward_reads), ignore.case = TRUE, perl = TRUE)
  }
  if (!is.character(sample_names) || anyNA(sample_names) || any(!nzchar(sample_names)) ||
      length(sample_names) != length(forward_reads) || anyDuplicated(sample_names)) {
    stop("sample_names must contain one unique name per read pair.", call. = FALSE)
  }

  for (i in seq_along(forward_reads)) validate_paired_fastq(forward_reads[i], reverse_reads[i])
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  filtered_dir <- file.path(output_dir, "filtered")
  dir.create(filtered_dir)
  file_names <- gsub("[^A-Za-z0-9_.-]", "_", sample_names)
  if (anyDuplicated(file_names)) {
    stop("Sample names produce duplicate file names.", call. = FALSE)
  }
  filtered_f <- file.path(filtered_dir, paste0(file_names, "_R1.fastq.gz"))
  filtered_r <- file.path(filtered_dir, paste0(file_names, "_R2.fastq.gz"))
  names(filtered_f) <- names(filtered_r) <- sample_names

  filtered_count <- dada2::filterAndTrim(
    forward_reads, filtered_f, reverse_reads, filtered_r,
    truncLen = as.integer(trunc_len), maxEE = max_ee, truncQ = 2L,
    rm.phix = TRUE, compress = TRUE, multithread = n_threads
  )
  if (any(filtered_count[, 2] == 0)) {
    stop("No reads passed filtering for: ", paste(sample_names[filtered_count[, 2] == 0], collapse = ", "),
         ". Inspect read quality and filtering parameters.", call. = FALSE)
  }
  error_f <- dada2::learnErrors(filtered_f, multithread = n_threads)
  error_r <- dada2::learnErrors(filtered_r, multithread = n_threads)
  derep_f <- lapply(filtered_f, dada2::derepFastq)
  derep_r <- lapply(filtered_r, dada2::derepFastq)
  names(derep_f) <- names(derep_r) <- sample_names
  denoised_f <- dada2::dada(derep_f, err = error_f, multithread = n_threads)
  denoised_r <- dada2::dada(derep_r, err = error_r, multithread = n_threads)
  merged <- dada2::mergePairs(
    denoised_f, derep_f, denoised_r, derep_r,
    minOverlap = 12L, maxMismatch = 0L
  )
  sequence_table <- dada2::makeSequenceTable(merged)
  sequence_table <- dada2::removeBimeraDenovo(
    sequence_table, method = "consensus", multithread = n_threads
  )

  count_reads <- function(x) sum(dada2::getUniques(x))
  tracking <- data.frame(
    sample = sample_names,
    input = filtered_count[, 1],
    filtered = filtered_count[, 2],
    denoised_f = vapply(denoised_f, count_reads, numeric(1)),
    denoised_r = vapply(denoised_r, count_reads, numeric(1)),
    merged = vapply(merged, count_reads, numeric(1)),
    nonchimera = rowSums(sequence_table),
    check.names = FALSE
  )

  result <- list(sequence_table = sequence_table, tracking = tracking)
  utils::write.csv(tracking, file.path(output_dir, "read_tracking.csv"),
                   row.names = FALSE)
  saveRDS(result, file.path(output_dir, "dada2_result.rds"))
  result
}
