
#' Convert a METAeDNA FASTA file to a sequence table
#'
#' \code{fasta_to_table()} reads a METAeDNA FASTA file, extracts
#' sequence-level metadata from `"demultiplex"`, `"stat"`, or `"denoise"`
#' JSON blocks and converts it into a wide abundance table. Demultiplexed reads
#' contribute a count of one and identical sequences are collapsed. Later
#' pipeline stages use their stored sample counts directly.
#'
#' @param fasta_file (character) path to the input FASTA file. The file must
#'   contain JSON metadata blocks in the header produced by earlier steps such
#'   as \code{demultiplex}, \code{collapse_sequences}, or \code{denoise}.
#'
#' @param n_threads (integer) number of parallel workers used for parsing JSON
#'   blocks. Records are divided into continuous batches to reduce parallel
#'   scheduling and data-combination overhead.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Reads the FASTA file using \code{Biostrings::readDNAStringSet()}.
#'   \item Extracts the JSON metadata block from each header. This block may
#'         originate from:
#'         \itemize{
#'           \item \code{"demultiplex"} (one count per assigned read),
#'           \item \code{"stat"} (output of the dereplication module), or
#'           \item \code{"denoise"} (output of the denoising module).
#'         }
#'   \item Parses per-experiment and per-sample counts and status labels.
#'         If \code{"denoise"} → \code{"by_sample"} is present, it is used;
#'         otherwise scalar sample counts in the \code{"stat"} block are used.
#'         A demultiplex record without either later block contributes one read
#'         to its stored experiment and sample.
#'   \item Constructs a long-format table with columns:
#'         \code{experiment}, \code{sample}, \code{count}, \code{status},
#'         \code{ID}, and \code{SEQ}.
#'   \item Computes total counts per sequence–experiment pair.
#'   \item Produces a wide-format abundance/status table using
#'         \code{tidyr::pivot_wider()}, where each sample contributes a
#'         \code{count_*} and \code{status_*} column.
#'   \item Assigns a unique sequence identifier of the form
#'         \code{SEQ_0000001}, \code{SEQ_0000002}, etc.
#' }
#'
#' @return
#' A wide-format \code{tibble} containing merged sequence metadata. The table
#' includes:
#' \itemize{
#'   \item \code{experiment}: experiment name,
#'   \item \code{ID}: original sequence identifier,
#'   \item \code{SEQ}: the nucleotide sequence as a character string,
#'   \item \code{SEQ_ID}: a newly assigned unified sequence identifier,
#'   \item \code{count_*}: per-sample abundance values,
#'   \item \code{status_*}: per-sample labels such as \code{"head"},
#'         \code{"internal"}, or \code{"singleton"},
#'   \item \code{total_count}: total abundance per sequence–experiment pair.
#' }
#'
#' @importFrom jsonlite fromJSON
#' @importFrom Biostrings readDNAStringSet
#' @importFrom dplyr mutate group_by summarise ungroup relocate
#' @importFrom tidyr pivot_wider
#'
#' @examples
#' \dontrun{
#' tbl <- fasta_to_table(
#'   "analysis_results/denoised_sequences.fasta",
#'   n_threads = 4
#' )
#' head(tbl)
#' }
#'
#' @export
fasta_to_table <- function(fasta_file, n_threads = 1) {
  if (length(fasta_file) != 1L || !file.exists(fasta_file)) {
    stop("fasta_file does not exist.", call. = FALSE)
  }
  if (length(n_threads) != 1L || is.na(n_threads) || n_threads < 1) {
    stop("n_threads must be a positive integer.", call. = FALSE)
  }

  dna_set <- Biostrings::readDNAStringSet(fasta_file)
  if (!length(dna_set)) {
    stop("The FASTA file does not contain any sequences.", call. = FALSE)
  }
  sequences <- as.character(dna_set)
  headers <- names(dna_set)
  sequence_ids <- sub("[;{].*$", "", headers)
  json_start <- regexpr("{", headers, fixed = TRUE)
  json_text <- rep(NA_character_, length(headers))
  has_json <- json_start > 0L
  json_text[has_json] <- substring(headers[has_json], json_start[has_json])
  if (!any(has_json)) {
    stop("No METAeDNA JSON metadata were found in the FASTA headers.",
         call. = FALSE)
  }
  if (any(!has_json)) {
    warning(sum(!has_json), " FASTA record(s) without JSON metadata were skipped.",
            call. = FALSE)
  }

  parse_record <- function(i) {
    js <- tryCatch(
      jsonlite::fromJSON(json_text[[i]], simplifyVector = FALSE),
      error = function(e) {
        stop("Invalid JSON metadata in FASTA record '", sequence_ids[[i]],
             "': ", conditionMessage(e), call. = FALSE)
      }
    )

    if (!is.null(js$denoise$by_sample)) {
      metadata <- js$denoise$by_sample
      rows <- lapply(names(metadata), function(experiment) {
        samples <- metadata[[experiment]]
        data.frame(
          experiment = experiment,
          sample = names(samples),
          count = vapply(samples, function(x) as.numeric(x$count), numeric(1)),
          status = vapply(samples, function(x) as.character(x$status),
                          character(1)),
          stringsAsFactors = FALSE
        )
      })
    } else if (!is.null(js$stat$exp)) {
      metadata <- js$stat$exp
      rows <- lapply(names(metadata), function(experiment) {
        samples <- metadata[[experiment]]$samples
        data.frame(
          experiment = experiment,
          sample = names(samples),
          count = as.numeric(unlist(samples, use.names = FALSE)),
          status = NA_character_,
          stringsAsFactors = FALSE
        )
      })
    } else if (!is.null(js$demultiplex)) {
      sample <- as.character(js$demultiplex$reaction_id)
      experiment <- as.character(js$demultiplex$exp)
      if (is.na(sample) || !nzchar(sample)) sample <- "no_sample"
      if (!length(experiment) || is.na(experiment) || !nzchar(experiment)) {
        experiment <- "no_experiment"
      }
      rows <- list(data.frame(
        experiment = experiment,
        sample = sample,
        count = 1,
        status = NA_character_,
        stringsAsFactors = FALSE
      ))
    } else {
      return(NULL)
    }

    result <- do.call(rbind, rows)
    if (is.null(result) || !nrow(result)) return(NULL)
    result$ID <- sequence_ids[[i]]
    result$SEQ <- unname(sequences[[i]])
    result
  }

  valid_index <- which(has_json)
  worker_count <- min(as.integer(n_threads), length(valid_index))
  chunk_id <- ceiling(seq_along(valid_index) * worker_count /
                        length(valid_index))
  index_chunks <- split(valid_index, chunk_id)
  parse_chunk <- function(index) {
    records <- lapply(index, parse_record)
    records <- records[!vapply(records, is.null, logical(1))]
    if (!length(records)) return(NULL)
    do.call(rbind, records)
  }

  if (worker_count == 1L) {
    parsed_chunks <- lapply(index_chunks, parse_chunk)
  } else {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = worker_count)
    parsed_chunks <- future.apply::future_lapply(
      index_chunks, parse_chunk, future.seed = TRUE
    )
  }
  parsed_chunks <- parsed_chunks[!vapply(parsed_chunks, is.null, logical(1))]
  if (!length(parsed_chunks)) {
    stop("No demultiplex, stat, or denoise sample metadata were found.",
         call. = FALSE)
  }
  all_seq_df <- do.call(rbind, parsed_chunks)
  rownames(all_seq_df) <- NULL

  first_id <- !duplicated(all_seq_df$SEQ)
  id_by_sequence <- stats::setNames(all_seq_df$ID[first_id],
                                    all_seq_df$SEQ[first_id])
  all_seq_df$ID <- unname(id_by_sequence[all_seq_df$SEQ])
  all_seq_df <- all_seq_df %>%
    group_by(SEQ, ID, experiment, sample) %>%
    summarise(
      count = sum(count),
      status = {
        value <- unique(stats::na.omit(status))
        if (length(value)) value[[1]] else NA_character_
      },
      .groups = "drop"
    ) %>%
    group_by(SEQ, experiment) %>%
    mutate(total_count = sum(count)) %>%
    ungroup()

  by_sample_df <- all_seq_df %>%
    tidyr::pivot_wider(
      id_cols = c(experiment, ID, SEQ, total_count),
      names_from = sample,
      values_from = c(count, status),
      values_fill = list(count = 0, status = NA_character_),
      names_sep = "_"
    )
  unique_sequences <- unique(by_sample_df$SEQ)
  by_sample_df$SEQ_ID <- sprintf(
    "SEQ_%07d", match(by_sample_df$SEQ, unique_sequences)
  )
  by_sample_df %>% relocate(experiment, ID, SEQ, total_count)
}

# Format Genbank fasta to dada2 style ---------------

#' Format GenBank FASTA records to the DADA2 reference format
#'
#' \code{format_to_dada()} converts a GenBank-style FASTA file into a DADA2-compatible
#' reference FASTA, where sequence names follow the standardized
#' \code{Kingdom;Phylum;Class;Order;Family;Genus;Species} taxonomy format.
#' Accession numbers are parsed from the original FASTA headers, taxonomic
#' identifiers are retrieved from a local SQL database, and complete taxonomy
#' strings are constructed for use in DADA2 reference workflows.
#'
#' @param fasta_file (character) path to the input FASTA file. The file
#'   should contain GenBank-style headers from which accession numbers can be
#'   extracted.
#'
#' @param output_fasta (character) path to the output FASTA file written
#'   in DADA2 reference format. Sequence headers are replaced by taxonomy
#'   strings with spaces converted to underscores.
#'
#' @param sql_path (character) path to a SQLite database containing NCBI
#'   taxonomy tables. This database must support accession → taxid mapping
#'   and taxid → full taxonomy queries used by \code{get_taxid()} and
#'   \code{get_taxonomy()}.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Reads all sequences using \code{Biostrings::readDNAStringSet()}.
#'   \item Extracts accession numbers from each GenBank-style header. Legacy
#'         \code{gi|} records are parsed individually from their embedded
#'         \code{"gb|ACC|"} field, so legacy and ordinary accession headers may
#'         be mixed in the same FASTA file.
#'   \item Queries a local SQL taxonomy database to obtain:
#'         \itemize{
#'           \item taxonomic identifiers using \code{get_taxid()},
#'           \item full taxonomy ranks using \code{get_taxonomy()}.
#'         }
#'   \item Constructs the DADA2-style taxonomy string:
#'         \code{kingdom;phylum;class;order;family;genus;species}.
#'   \item Replaces original sequence names with the constructed taxonomy
#'         strings (spaces are replaced by underscores).
#'   \item Writes the formatted sequences to \code{output_fasta} using
#'         \code{Biostrings::writeXStringSet()}.
#' }
#' Elapsed time is reported when the function exits.
#'
#' @return
#' Invisibly returns \code{NULL}. The output FASTA file contains the original
#' sequences with headers rewritten in DADA2 taxonomy format.
#'
#' @importFrom dplyr mutate as_tibble group_by reframe left_join rename_with ungroup relocate
#' @importFrom stringr str_extract str_remove str_detect str_replace_all
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#'
#' @examples
#' \dontrun{
#' # Use an accession-labelled FASTA and the local NCBI taxonomy database.
#' format_to_dada("reference.fasta", "dada2_reference.fasta", "ncbi.sqlite")
#' }#'
#' @export
format_to_dada <- function(fasta_file, output_fasta, sql_path){
  start_time <- Sys.time()
  on.exit(show_time(secs = as.numeric(Sys.time() - start_time, units = "secs"), mode = 1))
  raw_seq <- readDNAStringSet(fasta_file)

  seq_acc <- raw_seq %>% names() %>% str_remove(" .*")
  gi_rows <- str_detect(seq_acc, "^gi[|]")
  if(any(gi_rows)){
    parsed_gi <- seq_acc[gi_rows] %>%
      str_extract("gb[|][^|]+[|]") %>%
      str_remove("^gb[|]") %>%
      str_remove("[|]$")
    seq_acc[gi_rows] <- parsed_gi
  }
  taxonomy_id <- get_taxid(sql_path = sql_path, seq_acc)

  query_res <- get_taxonomy(sql_path = sql_path, taxids = taxonomy_id$taxid %>% unique, dereplicate = F)

  taxonomy_data <- left_join(taxonomy_id, query_res, by = join_by(taxid == query_taxid)) %>%
    mutate(dada_name = paste(kingdom, phylum, class, order, family, genus, species, sep = ";"))
  names(raw_seq) <- taxonomy_data$dada_name %>% str_replace_all(" ", "_")
  show_time("Formatting names...")
  writeXStringSet(raw_seq, output_fasta)
}
