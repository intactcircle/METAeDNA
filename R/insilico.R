#' Find primer binding sites on the target DNA sequence(s)
#'
#' This function identifies potential primer binding positions on a target DNA
#' sequence using approximate string matching. It supports forward and reverse
#' primers and allows a user-specified number of mismatches when searching for
#' primer–template alignments.
#'
#' @param target_sequence (character) path to a FASTA file containing the DNA
#'   sequence(s) to be searched. The file is read using
#'   \code{Biostrings::readDNAStringSet()}.
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
#' @details
#' The function:
#' \itemize{
#'   \item reads the target sequence(s) into a \code{DNAStringSet},
#'   \item optionally reverse-complements the primer,
#'   \item performs approximate pattern matching using
#'         \code{Biostrings::vmatchPattern()},
#'   \item returns a list-like object indicating all match positions for each
#'         input sequence.
#' }
#'
#' @return
#' A \code{MatchesList} object (from the Biostrings package) giving the locations
#' of primer binding sites within each target sequence. Each element contains a
#' \code{IRanges} object describing the positions of matches.
#'
#' @importFrom Biostrings readDNAStringSet DNAString reverseComplement vmatchPattern
#'
#' @examples
#' \dontrun{
#' # Search for a forward primer allowing up to 2 mismatches
#' find_primer_binding_sites(
#'   target_sequence = "chloroplast.fasta",
#'   primer = "AGCTTAGGCTAC",
#'   is_reverse = FALSE,
#'   max_mismatch = 2
#' )
#'
#' # Search for a reverse primer by providing its forward orientation
#' find_primer_binding_sites(
#'   target_sequence = "chloroplast.fasta",
#'   primer = "TCGATACCGGTA",
#'   is_reverse = TRUE,
#'   max_mismatch = 3
#' )
#' }
#'
#' @export
find_primer_binding_sites <- function(target_sequence, primer, is_reverse = FALSE, max_mismatch = 3) {
  # 转换目标序列为 DNAString 对象
  target_sequence <- readDNAStringSet(target_sequence)
  # 如果是反向引物，生成反向互补序列
  if (is_reverse) {
    primer <- reverseComplement(DNAString(primer))
  } else {
    primer <- DNAString(primer)
  }

  # 使用 vmatchPattern 查找匹配位置，允许一定的错配
  matches <- vmatchPattern(primer, target_sequence, max.mismatch = max_mismatch)
  # 返回结果
  return(matches)
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
#' This function performs an in-silico PCR search on a reference database by
#' locating primer-binding sites for forward and reverse primers, allowing
#' mismatches, optional positional constraints, and optional trimming of primers
#' from the resulting amplicons. The function supports degenerate bases and
#' allows users to embed constraint markers (`#`) inside primers to enforce
#' exact matching at specified positions. The search is based on approximate
#' pattern matching using \code{Biostrings::vmatchPattern()}.
#'
#' @param ref_db A \code{DNAStringSet} containing reference sequences to be
#'   searched.
#'
#' @param forward (character) Forward primer sequence in 5'→3' orientation.
#'   Degenerate IUPAC bases are allowed. If the sequence contains one or more
#'   `#` characters, these positions are treated as *strict* constraints:
#'   after approximate matching, the corresponding positions in the reference
#'   must match exactly.
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
#'   \item Optionally performs a reverse-strand search (if the object \code{try_reverse}
#'         is set in the calling environment; not user-exposed here).
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
#' The names of the returned sequences correspond to the reference IDs and
#' exclusion of the internal index used during processing.
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
#' @export
in_silico_base <- function(ref_db, forward, reverse, mis_match, trim, max_hits, ...) {
  forward_raw <- forward
  reverse_raw <- reverse
  constraint <- F
  #browser()
  # find restrictions
  if (str_detect(forward, "#")) {
    forward_rst <- (forward %>% gregexpr("#", .))[[1]] - 1
    constraint <- T
  }
  if (str_detect(reverse_raw, "#")) {
    reverse_rst <- (reverse_raw %>% gregexpr("#", .))[[1]]
    reverse_rst <- ((nchar(reverse_raw) - reverse_rst) + 1) %>% sort()
    constraint <- T
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
  # browser()
  if(nrow(match_rev_table) == 0 |nrow(match_fwd_table) == 0 ){
    ref_seqs <- DNAStringSet()
  }else{
    # generating a range table to extract sequences
    range_table <- full_join(match_rev_table, match_fwd_table,
                             relationship = "many-to-many", by = join_by(names)
    ) %>%
      filter(!is.na(rev_width), !is.na(fwd_width)) %>%
      mutate(
        temp_fwd_end = ifelse(fwd_end > rev_end, rev_end, fwd_end),
        temp_rev_end = ifelse(fwd_end > rev_end, fwd_end, rev_end),
        temp_fwd_start = ifelse(fwd_end > rev_end, rev_start, fwd_start),
        temp_rev_start = ifelse(fwd_end > rev_end, fwd_start, rev_start)
      ) %>%
      select(-fwd_start, -fwd_end, -rev_start, -rev_end) %>%
      dplyr::rename(
        fwd_start = temp_fwd_start,
        fwd_end = temp_fwd_end,
        rev_start = temp_rev_start,
        rev_end = temp_rev_end
      ) %>%
      filter(fwd_end < rev_start)

    # find constraints -------------------
    if (constraint) {
      strict_out <- c()
      temp_count <- 0
      ## forward constraint
      for (n in forward_rst) {
        pos <- nchar(forward) - (n - temp_count)
        subset_seqs <- subseq(ref_db[range_table$names], start = range_table$fwd_end - pos, end = range_table$fwd_end - pos) %>%
          as.character()
        check_chr <- forward %>% str_sub(start = n - temp_count, end = n - temp_count)
        #
        strict_out <- c(strict_out, subset_seqs[!subset_seqs %>% str_detect(convert_to_regex(check_chr))] %>% names())
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
        #
        strict_out <- c(strict_out, subset_seqs[!subset_seqs %>% str_detect(convert_to_regex(check_chr))] %>% names())
        temp_count <- temp_count + 1
      }
      range_table <- range_table %>% filter(!names %in% strict_out)
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
    names(ref_db) <- names(ref_db) %>% paste0("[REVERSE]",.)
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
    # browser()
    if(nrow(match_rev_table) == 0 |nrow(match_fwd_table) == 0 ){
      ref_seqs_rev <- DNAStringSet()
    }else{
      # generating a range table to extract sequences
      range_table <- full_join(match_rev_table, match_fwd_table,
                               relationship = "many-to-many", by = join_by(names)
      ) %>%
        filter(!is.na(rev_width), !is.na(fwd_width)) %>%
        mutate(
          temp_fwd_end = ifelse(fwd_end > rev_end, rev_end, fwd_end),
          temp_rev_end = ifelse(fwd_end > rev_end, fwd_end, rev_end),
          temp_fwd_start = ifelse(fwd_end > rev_end, rev_start, fwd_start),
          temp_rev_start = ifelse(fwd_end > rev_end, fwd_start, rev_start)
        ) %>%
        select(-fwd_start, -fwd_end, -rev_start, -rev_end) %>%
        dplyr::rename(
          fwd_start = temp_fwd_start,
          fwd_end = temp_fwd_end,
          rev_start = temp_rev_start,
          rev_end = temp_rev_end
        ) %>%
        filter(fwd_end < rev_start)

      # find constraints -------------------
      if (constraint) {
        strict_out <- c()
        temp_count <- 0
        ## forward constraint
        for (n in forward_rst) {
          pos <- nchar(forward) - (n - temp_count)
          subset_seqs <- subseq(ref_db[range_table$names], start = range_table$fwd_end - pos, end = range_table$fwd_end - pos) %>%
            as.character()
          check_chr <- forward %>% str_sub(start = n - temp_count, end = n - temp_count)
          #
          strict_out <- c(strict_out, subset_seqs[!subset_seqs %>% str_detect(convert_to_regex(check_chr))] %>% names())
          temp_count <- temp_count + 1
        }
        ## reverse constraint
        # CCGTCAATT#HC#TTY#AAR
        temp_count <- 0
        for (n in reverse_rst) {
          subset_seqs <- subseq(ref_db[range_table$names], start = range_table$rev_start + n - 1, end = range_table$rev_start + n - 1) %>%
            as.character()
          check_chr <- reverse %>%
            DNAString() %>%
            reverseComplement() %>%
            as.character() %>%
            str_sub(start = n - temp_count, end = n - temp_count)
          #
          strict_out <- c(strict_out, subset_seqs[!subset_seqs %>% str_detect(convert_to_regex(check_chr))] %>% names())
          temp_count <- temp_count + 1
        }
        range_table <- range_table %>% filter(!names %in% strict_out)
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
#' This function performs an in-silico PCR scan on a potentially large FASTA database
#' by iteratively loading the reference sequences in chunks, locating primer-binding
#' sites with mismatch tolerance, applying optional strict positional constraints
#' (via `#` markers), and extracting predicted amplicons. The computation is
#' delegated to \code{in_silico_base()}, but with chunked I/O to reduce memory usage.
#'
#' @param fasta_file (character) Path to a FASTA file containing reference sequences.
#'   The file may be large; sequences are loaded in chunks to reduce memory usage.
#'
#' @param forward (character) Forward primer in 5'→3' orientation. Degenerate bases
#'   are supported. Primer sequences may include `#` markers to indicate positional
#'   *exact-match constraints*. For example, \code{"CCGTCAATT#HC#TTY#AAR"} specifies
#'   that bases preceding `#` must match exactly in addition to general mismatch
#'   tolerance.
#'
#' @param reverse (character) Reverse primer in 5'→3' orientation (before
#'   reverse-complementing). Constraint markers (`#`) follow the same rules as in
#'   the forward primer and are evaluated after reverse-complementing.
#'
#' @param mis_match (integer >= 0) Maximum number of mismatches allowed during primer
#'   binding detection. Indels are not permitted.
#'
#' @param try_reverse (logical) If \code{TRUE}, an additional search is performed in
#'   which the entire reference sequence is reverse-complemented prior to matching.
#'   Returned amplicons from this mode are labeled with the tag \code{"[REVERSE]"}.
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
#' @param keep_genome_threshhold (integer >= 1) Length threshold defining a
#'   "genome-like" sequence. Any reference sequence with length >= this value
#'   will be excluded when \code{keep_genome = FALSE}.
#'
#' @param show_progress (logical) If \code{TRUE}, displays a text progress bar as
#'   chunks are processed.
#'
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @details
#' This function uses \code{\link[METAeDNA]{in_silico_base()}} as its core computational engine but
#' controls memory usage by loading only subsets of sequences at each iteration.
#'
#' The workflow is:
#'
#' \enumerate{
#'   \item Index the FASTA file using \code{\link[Biostrings]{fasta.index()}}.
#'   \item Optionally filter long (genomic) sequences for performance.
#'   \item Process the database in chunks of \code{chunk_size}.
#'   \item For each chunk, call \code{\link[METAeDNA]{in_silico_base()}} to:
#'     \itemize{
#'       \item locate approximate primer matches using \code{\link[Biostrings]{vmatchPattern()}},
#'       \item apply strict `#` constraints, if present,
#'       \item ensure correct forward→reverse orientation,
#'       \item extract amplicons (trimmed or untrimmed),
#'       \item optionally evaluate reverse-strand amplification.
#'     }
#'   \item Concatenate results from all chunks into a single \code{\link[Biostrings]{DNAStringSet()}}.
#' }
#'
#' This enables efficient simulation of PCR amplification against multi-million
#' sequence reference databases that cannot be fully loaded into memory.
#'
#' @return
#' A \code{DNAStringSet} containing all predicted amplicons across all processed
#' reference sequences. Names correspond to reference sequence identifiers, with
#' an optional \code{"[REVERSE]"} prefix when applicable.
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet fasta.index
#' @importFrom dplyr filter
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' \dontrun{
#' # Simple in-silico PCR simulation on a large database
#' res <- in_silico(
#'   fasta_file = "ref_database.fasta",
#'   forward = "CCGTCAATT#HC#TTY#AAR",
#'   reverse = "GACTACHVGGGTATCTAATCC",
#'   mis_match = 2,
#'   trim = TRUE,
#'   chunk_size = 5000,
#'   try_reverse = TRUE
#' )
#'
#' res
#' }
#' @export
in_silico <-
  function(fasta_file, forward, reverse, mis_match = 3, try_reverse = F, trim = T,chunk_size = 10000, max_hits = 5, keep_genome = F, keep_genome_threshhold = 10000,show_progress = FALSE, ...){
    fas_index <- fasta.index(fasta_file, seqtype="DNA")
    geno <- fas_index %>% filter(seqlength >= keep_genome_threshhold)
    if(nrow(geno) !=0 & keep_genome == F){
      warning("Genomic sequences presented. They are removed by default to speed up processing. To retain them, please set keep_genome = TRUE.")
    }
    if(!keep_genome){
      fas_index <- fas_index %>% filter(seqlength < keep_genome_threshhold)
    }
    total_rows <- fas_index %>% nrow()

    num_chunks <- ceiling(total_rows / chunk_size)
    if (show_progress) pb <- txtProgressBar(min = 0, max = num_chunks, style = 3)
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
        trim = trim,
        max_hits = max_hits
      )
      #browser()
      temp_res <- do.call(in_silico_base, args)
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
#' This function performs in-silico PCR for multiple primer pairs in parallel.
#' For each forward–reverse primer pair, it calls \code{in_silico()} and returns
#' a list of predicted amplicons. The computation is optionally parallelized
#' using the \pkg{future} and \pkg{future.apply} frameworks, and includes
#' real-time progress reporting via \pkg{progressr}.
#'
#' The function is designed for large FASTA databases, where each primer pair
#' must be evaluated independently and efficiently.
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
#'   labeled \code{"[REVERSE]"}.
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
#' @param workers (integer)
#'   Number of parallel workers.
#'   \code{workers = 1} → sequential execution.
#'   \code{workers > 1} → parallel execution via \code{future::multisession()}.
#'
#' @param ... Additional arguments forwarded to \code{\link[METAeDNA]{in_silico()}}.
#'
#' @return
#' A **named list** of \code{\link[Biostrings]{DNAStringSet()}} objects, one per primer pair.
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
#'   \item Use \code{\link[future.apply]{future_lapply()}} to loop over primer pairs.
#'   \item For each pair, call \code{\link[METAeDNA]{in_silico()}} to:
#'     \itemize{
#'       \item detect degenerate primer matches,
#'       \item apply positional \code{"#"} restrictions,
#'       \item identify valid amplicons,
#'       \item optionally include reverse-strand amplification.
#'     }
#'   \item Collect all results into a named list.
#' }
#'
#' The function uses \pkg{progressr} to report real-time progress when running
#' multiple primer pairs.
#'
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom progressr with_progress progressor handlers
#' @importFrom Biostrings DNAString DNAStringSet matchPattern readDNAStringSet reverseComplement vmatchPattern
#'
#' @examples
#' \dontrun{
#' res <- in_silico_multiprimer(
#'   fasta_file = "mydb.fasta",
#'   forwards   = c("CCGTCAATT#HC#TTY#AAR", "GACTCCTACGGGAGGC"),
#'   reverses   = c("GACTACHVGGGTATCTAATCC", "GTATTACCGCGGCTGCT"),
#'   mis_match  = 2,
#'   try_reverse = TRUE,
#'   workers    = 4
#' )
#'
#' names(res)
#' length(res[[1]])
#' }
#' @export
in_silico_multiprimer <- function(
    fasta_file,
    forwards,
    reverses,
    mis_match = 3,
    try_reverse = FALSE,
    trim = TRUE,
    chunk_size = 10000,
    max_hits = 5,
    workers = 1,
    ...) {
  # sanity check
  if (length(forwards) != length(reverses)) {
    stop("The forward and reverse primer numbers are inconsistent!")
  }

  # set parallel plan
  if (workers > 1) {
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::sequential)
  }

  # progress setup
  progressr::handlers(global = TRUE)
  progressr::handlers("progress")
  start_time <- format(Sys.time(), "%H:%M:%S")
  n_primers <- length(forwards)
  message("[",start_time  ,"] Running ", n_primers, " primer pairs on ", workers, " worker(s)...")

  # assign readable names
  primer_names <- sprintf("F_%s__R_%s", forwards, reverses)

  results <- NULL
  progressr::with_progress({
    p <- progressr::progressor(steps = n_primers)
    results <- future.apply::future_lapply(seq_len(n_primers), function(i) {
      p()
      in_silico(
        fasta_file = fasta_file,
        forward = forwards[i],
        reverse = reverses[i],
        mis_match = mis_match,
        try_reverse = try_reverse,
        trim = trim,
        chunk_size = chunk_size,
        max_hits = max_hits,
        ...
      )
    }, future.seed = T)
  })

  names(results) <- primer_names
  future::plan(future::sequential)
  return(results)
}



