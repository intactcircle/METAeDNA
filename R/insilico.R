#' Find primer binding sites in DNA sequences
#'
#' Search a FASTA file or DNA sequences for primer-binding positions. Reverse
#' primers are reverse-complemented before matching.
#'
#' @param target_sequence A FASTA file path, a DNA sequence or character vector
#'   of sequences, or a \code{DNAString}/\code{DNAStringSet} object. An existing
#'   file path is read as FASTA; character values otherwise represent sequences.
#'
#' @param primer (character) primer sequence in 5'→3' orientation. If
#'   \code{is_reverse = TRUE}, the function automatically computes the reverse
#'   complement before matching.
#'
#' @param is_reverse (logical) whether the provided primer should be treated as a
#'   reverse primer. When \code{TRUE}, the reverse complement of the primer is
#'   matched to the target sequence.
#'
#' @param max_mismatch (integer) maximum number of mismatched bases allowed in
#'   the approximate match, passed to \code{vmatchPattern()} via
#'   \code{max.mismatch}. Higher values broaden the search and may yield multiple
#'   candidate binding regions.
#'
#' @return A data frame with one row per binding site: `sequence_name`, `start`,
#'   `end`, and `strand`. Positions are 1-based and include both primer ends.
#'   Unnamed sequences are labelled `seq_1`, `seq_2`, and so on. An empty data
#'   frame with the same columns is returned when no sites match.
#'
#' @importFrom Biostrings readDNAStringSet DNAString DNAStringSet reverseComplement vmatchPattern
#'
#' @examples
#' find_primer_binding_sites("ACGTACTGATCGATGCTAGCTAGC",
#'                           primer = "ACGTACTG", max_mismatch = 0)
#' find_primer_binding_sites(c(a = "ACGTACTG", b = "TTTACGTACTG"),
#'                           primer = "ACGTACTG", max_mismatch = 0)
#' fasta <- tempfile(fileext = ".fasta")
#' writeLines(c(">reference", "ACGTACTGATCGATGCTAGCTAGC"), fasta)
#' find_primer_binding_sites(fasta, primer = "ACGTACTG", max_mismatch = 0)
#' unlink(fasta)
#' @export
find_primer_binding_sites <- function(target_sequence, primer, is_reverse = FALSE, max_mismatch = 3) {
  if (is.character(target_sequence) && length(target_sequence) == 1L &&
      file.exists(target_sequence)) {
    target_sequence <- readDNAStringSet(target_sequence)
  } else if (!inherits(target_sequence, "DNAStringSet")) {
    target_sequence <- DNAStringSet(target_sequence)
  }

  primer <- DNAString(primer)
  if (is_reverse) primer <- reverseComplement(primer)
  matches <- vmatchPattern(primer, target_sequence, max.mismatch = max_mismatch)
  positions <- as.data.frame(matches)
  if (!nrow(positions)) {
    return(data.frame(sequence_name = character(), start = integer(),
                      end = integer(), strand = character()))
  }

  sequence_names <- names(target_sequence)
  if (is.null(sequence_names)) sequence_names <- rep("", length(target_sequence))
  unnamed <- is.na(sequence_names) | !nzchar(sequence_names)
  sequence_names[unnamed] <- paste0("seq_", which(unnamed))
  data.frame(sequence_name = sequence_names[positions$group],
             start = positions$start, end = positions$end,
             strand = if (is_reverse) "-" else "+")
}
# create regex for DNA coding  ---------------
convert_to_regex <- function(dna_seq) {
  # 定义简并碱基
  degenerate_bases <- c(
    "A" = "A", "T" = "T", "C" = "C", "G" = "G",
    "R" = "[AG]", "Y" = "[CT]", "S" = "[GC]", "W" = "[AT]",
    "K" = "[GT]", "M" = "[AC]", "B" = "[CGT]", "D" = "[AGT]",
    "H" = "[ACT]", "V" = "[ACG]", "N" = "[ACGT]"
  )

  # 替换 DNA 序列中的简并碱基
  regex_seq <- degenerate_bases[dna_seq]

  return(regex_seq)
}

# in silico pcr base function ----------------------------------
#' In-silico PCR: identify primer-binding sites and extract predicted amplicons
#'
#' \code{in_silico_base()} searches a reference database by
#' locating primer-binding sites for forward and reverse primers, allowing
#' mismatches, optional positional constraints, and optional trimming of primers
#' from the resulting amplicons. The function supports degenerate bases and
#' accepts constraint markers (`#`) inside primers to require
#' exact matching at specified positions. The search is based on approximate
#' pattern matching using \code{Biostrings::vmatchPattern()}.
#'
#' @param ref_db A \code{DNAStringSet} containing reference sequences to be
#'   searched.
#'
#' @param forward (character) Forward primer sequence in 5'→3' orientation.
#'   Degenerate IUPAC bases are allowed. If the sequence contains one or more
#'   `#` characters, the primer base immediately before each marker must match
#'   exactly. The marker does not constrain the following base and is removed
#'   before primer matching.
#'
#' @param reverse (character) Reverse primer sequence in 5'→3' orientation before
#'   reverse-complementing. The function automatically computes the reverse
#'   complement for matching. Constraint markers (`#`) are interpreted relative
#'   to the original forward orientation before complementing.
#'
#' @param mis_match (integer >= 0) Maximum number of mismatches allowed when
#'   locating primer-binding sites using \code{vmatchPattern(max.mismatch = ...)}.
#'   Insertions/deletions are not allowed (Hamming-distance matching).
#'
#' @param trim (logical) If \code{TRUE}, the returned sequences contain only the
#'   internal amplicon region (primers removed). If \code{FALSE}, both primers
#'   are retained in the extracted sequence.
#'
#' @param max_hits (integer >= 1) Maximum number of primer-binding locations to
#'   retain per reference sequence for both forward and reverse primers. Only
#'   the first \code{max_hits} matches (ordered by genomic location) are used in
#'   constructing possible amplicons.
#' @param try_reverse (logical) also search the reverse complement of each
#'   reference sequence.
#' @param mark_reverse (logical) If \code{TRUE}, append the tag
#'   \verb{[REVERSE]} to names of amplicons produced by the reverse-complement
#'   search.
#'
#' @param ... Additional arguments (currently ignored; included for forward
#'   compatibility).
#'
#' @details
#' The function:
#'
#' \enumerate{
#'   \item Removes \code{"#"} from the primer sequences while storing their
#'         original constrained positions.
#'   \item Matches forward and reverse-complement primers using
#'         \code{vmatchPattern()} with \code{max.mismatch = mis_match}.
#'   \item Enumerates all valid forward–reverse match pairs where
#'         \code{fwd\_end < rev\_start}.
#'   \item Applies strict positional constraints at `#` locations, ensuring
#'         that the reference sequence base exactly matches the primer base.
#'   \item Optionally trims primers from the extracted region.
#'   \item Optionally performs a reverse-strand search when
#'         \code{try_reverse = TRUE}.
#' }
#'
#' This procedure simulates PCR amplification logic: both primers must bind in
#' the correct orientation, on the same sequence, and in forward → reverse order.
#'
#' @return
#' A \code{DNAStringSet} containing all predicted amplicons that satisfy:
#'
#' \itemize{
#'   \item allowed mismatch thresholds,
#'   \item correct forward–reverse ordering,
#'   \item optional strict `#` constraints,
#'   \item optional primer trimming.
#' }
#'
#' Returned sequence names retain the reference IDs without the processing index.
#'
#' @importFrom Biostrings DNAString DNAStringSet reverseComplement vmatchPattern subseq width
#' @importFrom stringr str_detect str_remove_all str_remove str_sub
#' @importFrom dplyr group_by slice_head ungroup full_join rename mutate filter left_join
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#'
#' ref <- DNAStringSet(c(
#'   seq1 = "ACGTACGTACGTACGTACGTACGT",
#'   seq2 = "TTTTACGTACGTTTTACGTACGTAA"
#' ))
#'
#' # Forward primer with a strict position (#)
#' forward <- "ACGTA#GT"
#'
#' # Reverse primer (before reverse complement)
#' reverse <- "ACGTACGT"
#'
#' in_silico_base(
#'   ref_db = ref,
#'   forward = forward,
#'   reverse = reverse,
#'   mis_match = 2,
#'   trim = TRUE,
#'   max_hits = 5
#' )
#' }
#'
#' @noRd
in_silico_base <- function(ref_db, forward, reverse, mis_match = 3, trim = TRUE,
                           max_hits = 5, try_reverse = FALSE,
                           mark_reverse = TRUE, ...) {
  if (grepl("^#|##", forward) || grepl("^#|##", reverse)) {
    stop("Each # must follow a primer base.", call. = FALSE)
  }
  reverse_raw <- reverse
  constraint <- FALSE
  forward_rst <- integer()
  reverse_rst <- integer()
  # find restrictions
  if (str_detect(forward, "#")) {
    forward_rst <- (forward %>% gregexpr("#", .))[[1]] - 1
    constraint <- TRUE
  }
  if (str_detect(reverse_raw, "#")) {
    reverse_rst <- (reverse_raw %>% gregexpr("#", .))[[1]]
    reverse_rst <- ((nchar(reverse_raw) - reverse_rst) + 1) %>% sort()
    constraint <- TRUE
  }

  # remove restrictions
  forward <- forward %>% str_remove_all("#")
  reverse <- reverse %>% str_remove_all("#")
  # read DNA string

  names(ref_db) <- names(ref_db) %>% paste0(., "[IDIND", 1:length(ref_db), "]")
  # string matching
  match_fwd <- vmatchPattern(forward, ref_db, max.mismatch = mis_match, fixed = FALSE)
  match_rev <- vmatchPattern(reverseComplement(reverse %>% DNAString()) %>% as.character(), ref_db, max.mismatch = mis_match, fixed = FALSE)
  # convert to table
  match_fwd_table <- unlist(match_fwd, recursive = TRUE, use.names = TRUE) %>%
    data.frame() %>%
    dplyr::rename(fwd_start = start, fwd_end = end, fwd_width = width) %>%
    group_by(names) %>%
    slice_head(n = max_hits) %>% ungroup
  match_rev_table <- unlist(match_rev, recursive = TRUE, use.names = TRUE) %>%
    data.frame() %>%
    dplyr::rename(rev_start = start, rev_end = end, rev_width = width) %>%
    group_by(names) %>%
    slice_head(n = max_hits)%>% ungroup
  if(nrow(match_rev_table) == 0 |nrow(match_fwd_table) == 0 ){
    ref_seqs <- DNAStringSet()
  }else{
    # generating a range table to extract sequences
    range_table <- full_join(match_rev_table, match_fwd_table,
                             relationship = "many-to-many", by = join_by(names)
    ) %>%
      filter(!is.na(rev_width), !is.na(fwd_width)) %>%
      filter(fwd_end < rev_start)

    # find constraints -------------------
    if (constraint) {
      strict_keep <- rep(TRUE, nrow(range_table))
      temp_count <- 0
      ## forward constraint
      for (n in forward_rst) {
        pos <- nchar(forward) - (n - temp_count)
        subset_seqs <- subseq(ref_db[range_table$names], start = range_table$fwd_end - pos, end = range_table$fwd_end - pos) %>%
          as.character()
        check_chr <- forward %>% str_sub(start = n - temp_count, end = n - temp_count)
        strict_keep <- strict_keep & str_detect(
          subset_seqs, convert_to_regex(check_chr)
        )
        temp_count <- temp_count + 1
      }
      ## reverse constraint
      # CCGTCAATT#HC#TTY#AAR
      temp_count <- 0
      for (n in reverse_rst) {
        subset_seqs <- subseq(ref_db[range_table$names], start = range_table$rev_start + n - 1 - temp_count, end = range_table$rev_start + n - temp_count -1) %>%
          as.character()
        check_chr <- reverse %>%
          DNAString() %>%
          reverseComplement() %>%
          as.character() %>%
          str_sub(start = n - temp_count, end = n - temp_count)
        strict_keep <- strict_keep & str_detect(
          subset_seqs, convert_to_regex(check_chr)
        )
        temp_count <- temp_count + 1
      }
      range_table <- range_table[strict_keep, , drop = FALSE]
    }
    #
    if (trim) {
      # remove primer
      ref_seqs <- subseq(ref_db[range_table$names], start = range_table$fwd_end + 1, end = range_table$rev_start - 1)
    } else {
      # keep primer
      range_table <- range_table %>%
        left_join(data.frame(names = names(ref_db), width = width(ref_db)), by = join_by(names)) %>%
        mutate(
          fwd_start = case_when(fwd_start <= 0 ~ 1,
                                .default = fwd_start
          ),
          rev_end = case_when(rev_end > width ~ width,
                              .default = rev_end
          )
        )
      ref_seqs <- subseq(ref_db[range_table$names], start = range_table$fwd_start, end = range_table$rev_end)
    }
  }

  # reverse match
  if(try_reverse){
    ref_db <- ref_db %>% reverseComplement()
    if (mark_reverse) names(ref_db) <- paste0(names(ref_db), "[REVERSE]")
    # string matching
    match_fwd <- vmatchPattern(forward, ref_db, max.mismatch = mis_match, fixed = FALSE)
    match_rev <- vmatchPattern(reverseComplement(reverse %>% DNAString()) %>% as.character(), ref_db, max.mismatch = mis_match, fixed = FALSE)
    # convert to table
    match_fwd_table <- unlist(match_fwd, recursive = TRUE, use.names = TRUE) %>%
      data.frame() %>%
      dplyr::rename(fwd_start = start, fwd_end = end, fwd_width = width)%>%
      group_by(names) %>%
      slice_head(n = max_hits)%>% ungroup
    match_rev_table <- unlist(match_rev, recursive = TRUE, use.names = TRUE) %>%
      data.frame() %>%
      dplyr::rename(rev_start = start, rev_end = end, rev_width = width)%>%
      group_by(names) %>%
      slice_head(n = max_hits)%>% ungroup
    if(nrow(match_rev_table) == 0 |nrow(match_fwd_table) == 0 ){
      ref_seqs_rev <- DNAStringSet()
    }else{
      # generating a range table to extract sequences
      range_table <- full_join(match_rev_table, match_fwd_table,
                               relationship = "many-to-many", by = join_by(names)
      ) %>%
        filter(!is.na(rev_width), !is.na(fwd_width)) %>%
        filter(fwd_end < rev_start)

      # find constraints -------------------
      if (constraint) {
        strict_keep <- rep(TRUE, nrow(range_table))
        temp_count <- 0
        ## forward constraint
        for (n in forward_rst) {
          pos <- nchar(forward) - (n - temp_count)
          subset_seqs <- subseq(ref_db[range_table$names], start = range_table$fwd_end - pos, end = range_table$fwd_end - pos) %>%
          as.character()
          check_chr <- forward %>% str_sub(start = n - temp_count, end = n - temp_count)
          strict_keep <- strict_keep & str_detect(
            subset_seqs, convert_to_regex(check_chr)
          )
          temp_count <- temp_count + 1
        }
        ## reverse constraint
        # CCGTCAATT#HC#TTY#AAR
        temp_count <- 0
        for (n in reverse_rst) {
          subset_seqs <- subseq(ref_db[range_table$names], start = range_table$rev_start + n - 1 - temp_count, end = range_table$rev_start + n - 1 - temp_count) %>%
            as.character()
          check_chr <- reverse %>%
            DNAString() %>%
            reverseComplement() %>%
            as.character() %>%
            str_sub(start = n - temp_count, end = n - temp_count)
          strict_keep <- strict_keep & str_detect(
            subset_seqs, convert_to_regex(check_chr)
          )
          temp_count <- temp_count + 1
        }
        range_table <- range_table[strict_keep, , drop = FALSE]
      }
      #
      if (trim) {
        # remove primer
        ref_seqs_rev <- subseq(ref_db[range_table$names], start = range_table$fwd_end + 1, end = range_table$rev_start - 1)
      } else {
        # keep primer
        range_table <- range_table %>%
          left_join(data.frame(names = names(ref_db), width = width(ref_db)), by = join_by(names)) %>%
          mutate(
            fwd_start = case_when(fwd_start <= 0 ~ 1,
                                  .default = fwd_start
            ),
            rev_end = case_when(rev_end > width ~ width,
                                .default = rev_end
            )
          )
        ref_seqs_rev <- subseq(ref_db[range_table$names], start = range_table$fwd_start, end = range_table$rev_end)
      }
      ref_seqs <- c(ref_seqs, ref_seqs_rev)
    }


  }
  #

  # remove ID indexes
  names(ref_seqs) <- names(ref_seqs) %>% str_remove("\\[IDIND[0-9]{1,}\\]")
  return(ref_seqs)
}

#' In-silico PCR Across Large FASTA Databases
#'
#' \code{in_silico()} scans a reference FASTA database
#' by iteratively loading the reference sequences in chunks, locating primer-binding
#' sites with mismatch tolerance, applying optional strict positional constraints
#' (via `#` markers), and extracting predicted amplicons. The computation is
#' performed with chunked I/O to reduce memory usage. One primer pair returns a
#' \code{DNAStringSet}; vectors of primer pairs return a named list and may be
#' processed in parallel.
#'
#' @param fasta_file (character) Path to a FASTA file containing reference sequences.
#'   The file may be large; sequences are loaded in chunks to reduce memory usage.
#'
#' @param forward (character vector) Forward primer or primers in 5'→3' orientation. Degenerate bases
#'   are supported. Primer sequences may include `#` markers to indicate positional
#'   *exact-match constraints*. For example, \code{"CCGTCAATT#HC#TTY#AAR"} specifies
#'   that the base immediately before each `#` must match exactly in addition to
#'   the general mismatch tolerance. The `#` marks the preceding base, not the
#'   following base, and is removed before primer matching.
#'
#' @param reverse (character vector) Reverse primer or primers in 5'→3'
#'   orientation, paired element-by-element with \code{forward}. Constraint
#'   markers (`#`) also constrain the base immediately before each marker in the
#'   supplied 5'→3' primer and are mapped to the corresponding position after
#'   reverse-complementing.
#'
#' @param mis_match (integer >= 0) Maximum number of mismatches allowed during primer
#'   binding detection. Indels are not permitted.
#'
#' @param try_reverse (logical) If \code{TRUE}, an additional search is performed in
#'   which the entire reference sequence is reverse-complemented prior to matching.
#' @param mark_reverse (logical) When \code{try_reverse = TRUE}, append
#'   \verb{[REVERSE]} to the end of names of reverse-search amplicons. Set it
#'   to \code{FALSE} to keep the original reference names. Default is
#'   \code{TRUE}.
#'
#' @param trim (logical) If \code{TRUE}, primers are removed from returned sequences;
#'   otherwise, primer bases are retained.
#'
#' @param chunk_size (integer >= 1) Number of sequences to load per iteration.
#'   Chunk-based processing reduces memory requirements for large reference
#'   databases.
#'
#' @param max_hits (integer >= 1) Maximum number of forward- or reverse-primer
#'   matches retained per sequence when building candidate amplicons.
#'
#' @param keep_genome (logical) If \code{FALSE} (default), very long sequences—
#'   typically genomic contigs—are removed before searching to greatly reduce
#'   computation time.
#'
#' @param keep_genome_threshold (integer >= 1) Length threshold defining a
#'   "genome-like" sequence. Any reference sequence with length >= this value
#'   will be excluded when \code{keep_genome = FALSE}.
#'
#' @param show_progress (logical) If \code{TRUE}, displays a text progress bar as
#'   chunks are processed.
#' @param n_threads Number of parallel R workers used when multiple primer pairs
#'   are supplied. A value of 1 processes the pairs sequentially.
#'
#' @param ... Reserved for future extensions; currently ignored.
#' @param min_length,max_length Optional minimum and maximum returned amplicon
#'   lengths in bases, inclusive. `NULL` leaves the corresponding bound unset.
#'   Length is measured after primer trimming when `trim = TRUE`, and includes
#'   both primers when `trim = FALSE`. These limits do not filter reference
#'   sequence lengths.
#'
#' @details
#' The function controls memory usage by loading only subsets of sequences at
#' each iteration.
#'
#' The workflow is:
#'
#' \enumerate{
#'   \item Index the FASTA file using \code{Biostrings::fasta.index()}.
#'   \item Optionally filter long (genomic) sequences for performance.
#'   \item Process the database in chunks of \code{chunk_size}.
#'   \item For each chunk:
#'     \itemize{
#'       \item locate approximate primer matches using \code{Biostrings::vmatchPattern()},
#'       \item apply strict `#` constraints, if present,
#'       \item require forward-to-reverse primer order,
#'       \item extract amplicons (trimmed or untrimmed),
#'       \item optionally evaluate reverse-strand amplification.
#'     }
#'   \item Concatenate results into a \code{Biostrings::DNAStringSet} object.
#' }
#'
#' The reference database is read in chunks rather than loaded in full.
#'
#' @return
#' For one primer pair, a \code{DNAStringSet} containing the predicted
#' amplicons. For multiple primer pairs, a named list of \code{DNAStringSet}
#' objects. Names correspond to reference identifiers, with an optional
#' \verb{[REVERSE]} suffix when applicable.
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet fasta.index width
#' @importFrom dplyr filter
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' \dontrun{
#' # Simple in-silico PCR simulation on a large database
#' res <- in_silico(
#'   fasta_file = "reference_databases/reference_sequences.fasta",
#'   forward = "CCGTCAATT#HC#TTY#AAR",
#'   reverse = "GACTACHVGGGTATCTAATCC",
#'   mis_match = 2,
#'   trim = TRUE,
#'   chunk_size = 5000,
#'   try_reverse = TRUE,
#'   mark_reverse = TRUE,
#'   min_length = 50,
#'   max_length = 500
#' )
#'
#' res
#'
#' # Multiple primer pairs
#' multi <- in_silico(
#'   fasta_file = "reference_databases/reference_sequences.fasta",
#'   forward = c("CCGTCAATT#HC#TTY#AAR", "GACTCCTACGGGAGGC"),
#'   reverse = c("GACTACHVGGGTATCTAATCC", "GTATTACCGCGGCTGCT"),
#'   n_threads = 2
#' )
#' }
#' @export
in_silico <-
  function(fasta_file = NULL, forward, reverse, mis_match = 3,
           try_reverse = FALSE, trim = TRUE,
           chunk_size = 10000, max_hits = 5, keep_genome = FALSE,
           keep_genome_threshold = 10000,
           show_progress = FALSE, n_threads = 1L,
           mark_reverse = TRUE, min_length = NULL, max_length = NULL, ...){
    if (is.null(fasta_file)) stop("fasta_file must be provided.", call. = FALSE)
    if (!is.null(min_length)) .check_number(min_length, "min_length", 1, .Machine$integer.max, TRUE)
    if (!is.null(max_length)) .check_number(max_length, "max_length", 1, .Machine$integer.max, TRUE)
    if (!is.null(min_length) && !is.null(max_length) && min_length > max_length) {
      stop("min_length must not exceed max_length.", call. = FALSE)
    }
    if (length(try_reverse) != 1L || is.na(try_reverse) || !is.logical(try_reverse)) {
      stop("try_reverse must be TRUE or FALSE.", call. = FALSE)
    }
    if (length(mark_reverse) != 1L || is.na(mark_reverse) || !is.logical(mark_reverse)) {
      stop("mark_reverse must be TRUE or FALSE.", call. = FALSE)
    }
    if (length(forward) != length(reverse)) {
      stop("forward and reverse must contain the same number of primers.", call. = FALSE)
    }
    if (length(forward) > 1L) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      if (n_threads > 1L) {
        future::plan(future::multisession, workers = n_threads)
      } else {
        future::plan(future::sequential)
      }
      results <- future.apply::future_lapply(seq_along(forward), function(i) {
        in_silico(
          fasta_file = fasta_file,
          forward = forward[i],
          reverse = reverse[i],
          mis_match = mis_match,
          try_reverse = try_reverse,
          mark_reverse = mark_reverse,
          trim = trim,
          chunk_size = chunk_size,
          max_hits = max_hits,
          keep_genome = keep_genome,
          keep_genome_threshold = keep_genome_threshold,
          show_progress = FALSE,
          n_threads = 1L,
          min_length = min_length,
          max_length = max_length,
          ...
        )
      }, future.seed = TRUE)
      names(results) <- sprintf("F_%s__R_%s", forward, reverse)
      return(results)
    }
    fas_index <- fasta.index(fasta_file, seqtype="DNA")
    geno <- fas_index %>% filter(seqlength >= keep_genome_threshold)
    if(nrow(geno) !=0 & keep_genome == FALSE){
      warning(
        "Genome-length sequences were removed; set keep_genome = TRUE to retain them.",
        call. = FALSE
      )
    }
    if(!keep_genome){
      fas_index <- fas_index %>% filter(seqlength < keep_genome_threshold)
    }
    total_rows <- fas_index %>% nrow()
    if (total_rows == 0L) return(DNAStringSet())

    num_chunks <- ceiling(total_rows / chunk_size)
    if (show_progress) pb <- txtProgressBar(min = 0, max = num_chunks, width = 40, style = 3)
    insilico_result <- DNAStringSet()
    for (i in seq_len(num_chunks)) {
      start <- (i - 1) * chunk_size + 1
      end <- min(i * chunk_size, total_rows)

      ref_db <- fas_index[start:end, ] %>% readDNAStringSet()
      args <- list(
        ref_db = ref_db,
        forward = forward,
        reverse = reverse,
        mis_match = mis_match,
        try_reverse = try_reverse,
        mark_reverse = mark_reverse,
        trim = trim,
        max_hits = max_hits
      )
      temp_res <- do.call(in_silico_base, args)
      lengths <- width(temp_res)
      keep <- rep(TRUE, length(temp_res))
      if (!is.null(min_length)) keep <- keep & lengths >= min_length
      if (!is.null(max_length)) keep <- keep & lengths <= max_length
      temp_res <- temp_res[keep]
      insilico_result <- c(insilico_result, temp_res)

      # 处理或保存 chunk

      if (show_progress) setTxtProgressBar(pb, i)
    }
    if (show_progress)close(pb)
    #
    return(insilico_result)
}


#' Parallel In-silico PCR for Multiple Primer Pairs
#'
#' \code{in_silico_multiprimer()} performs in-silico PCR for multiple primer pairs.
#' For each forward–reverse primer pair, it calls \code{in_silico()} and returns
#' a list of predicted amplicons. The computation is optionally parallelized
#' using the \pkg{future} and \pkg{future.apply} frameworks.
#'
#' Each primer pair is evaluated independently against the same database.
#'
#' @param fasta_file (character)
#'   Path to the reference FASTA file used for in-silico PCR.
#'   The file may be large; \code{in_silico()} internally loads sequences
#'   in chunks via \code{fasta.index()}.
#'
#' @param forwards (character vector)
#'   Forward primers in 5'→3'. Degenerate IUPAC bases are supported.
#'   Primer sequences may contain \code{"#"} markers to indicate exact-match
#'   positional constraints (see \code{in_silico()} for details).
#'
#' @param reverses (character vector)
#'   Reverse primers corresponding one-to-one with \code{forwards}.
#'   Constraint markers are supported and are applied after reverse-
#'   complementing the primer.
#'
#' @param mis_match (integer)
#'   Maximum number of allowed mismatches when locating primer binding sites.
#'
#' @param try_reverse (logical)
#'   If \code{TRUE}, \code{in_silico()} additionally tests amplification on
#'   the reverse complement of every reference sequence, returning amplicons
#'   from the opposite orientation.
#' @param mark_reverse (logical) If \code{TRUE}, append \verb{[REVERSE]} to
#'   names of reverse-search amplicons.
#'
#' @param trim (logical)
#'   If \code{TRUE}, primers are removed from returned amplicons.
#'
#' @param chunk_size (integer >= 1)
#'   Number of reference sequences processed per chunk. Passed directly to
#'   \code{in_silico()}. Helps reduce memory usage for large databases.
#'
#' @param max_hits (integer >= 1)
#'   Maximum number of primer-binding hits retained per reference sequence.
#'
#' @param n_threads (integer)
#'   Number of parallel workers.
#'   \code{n_threads = 1} → sequential execution.
#'   \code{n_threads > 1} → parallel execution via \code{future::multisession()}.
#'
#' @param ... Additional arguments forwarded to \code{\link{in_silico}}.
#'
#' @return
#' A named list of \code{Biostrings::DNAStringSet} objects, one per primer pair.
#' List names follow the pattern:
#'
#' \preformatted{
#'   "F_<forward>__R_<reverse>"
#' }
#'
#' @details
#' The workflow is:
#'
#' \enumerate{
#'   \item Verify that the number of forward and reverse primers matches.
#'   \item Initialize a parallel backend (or run sequentially).
#'   \item Use \code{future.apply::future_lapply()} to loop over primer pairs.
#'   \item For each pair, call \code{\link{in_silico}} to:
#'     \itemize{
#'       \item detect degenerate primer matches,
#'       \item apply positional \code{"#"} restrictions,
#'       \item identify valid amplicons,
#'       \item optionally include reverse-strand amplification.
#'     }
#'   \item Collect all results into a named list.
#' }
#'
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom Biostrings DNAString DNAStringSet matchPattern readDNAStringSet reverseComplement vmatchPattern
#'
#' @examples
#' \dontrun{
#' res <- in_silico_multiprimer(
#'   fasta_file = "reference_databases/reference_sequences.fasta",
#'   forwards   = c("CCGTCAATT#HC#TTY#AAR", "GACTCCTACGGGAGGC"),
#'   reverses   = c("GACTACHVGGGTATCTAATCC", "GTATTACCGCGGCTGCT"),
#'   mis_match  = 2,
#'   try_reverse = TRUE,
#'   n_threads    = 4
#' )
#'
#' names(res)
#' length(res[[1]])
#' }
#' @noRd
in_silico_multiprimer <- function(
    fasta_file,
    forwards,
    reverses,
    mis_match = 3,
    try_reverse = FALSE,
    mark_reverse = TRUE,
    trim = TRUE,
    chunk_size = 10000,
    max_hits = 5,
    n_threads = 1,
    ...) {
  # sanity check
  if (length(forwards) != length(reverses)) {
    stop("The forward and reverse primer numbers are inconsistent!")
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (n_threads > 1) {
    future::plan(future::multisession, workers = n_threads)
  } else {
    future::plan(future::sequential)
  }

  n_primers <- length(forwards)
  primer_names <- sprintf("F_%s__R_%s", forwards, reverses)

  results <- future.apply::future_lapply(seq_len(n_primers), function(i) {
    in_silico(
      fasta_file = fasta_file,
      forward = forwards[i],
      reverse = reverses[i],
      mis_match = mis_match,
      try_reverse = try_reverse,
      mark_reverse = mark_reverse,
      trim = trim,
      chunk_size = chunk_size,
      max_hits = max_hits,
      ...
    )
  }, future.seed = TRUE)

  names(results) <- primer_names
  results
}
