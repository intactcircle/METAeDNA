#' Infer Last Common Ancestor from an annotated BLAST table
#'
#' Calculate rank-wise LCA assignments from a BLAST data frame that has already
#' been annotated by \code{get_taxonomy()}. Input is an annotated data frame;
#' files and database paths are not accepted.
#'
#' @param df A data frame containing BLAST fields and taxonomy columns.
#' @param seqs_ID Name of the query-sequence column.
#' @param tax_levels Taxonomy columns ordered from the highest to the lowest
#'   rank.
#' @param taxid_column Name of the NCBI TaxID column.
#' @param accession_column Name of the subject accession column.
#' @param identity_column Name of the percentage-identity column.
#' @param coverage_column Name of the query-coverage column.
#' @param identity_threshold Optional named numeric vector giving the minimum
#'   percentage identity required at each rank, for example
#'   `c(species = 100, genus = 98, family = 95)`. The default `NULL` applies no
#'   identity threshold.
#'
#' @return A tibble with one LCA assignment per query sequence. If the input
#'   contains COL calibration columns, the result also reports whether all
#'   retained best hits were resolved and their accepted COL taxon IDs.
#'
#' @details
#' Use \code{get_taxonomy(blast_result = ...)} before this function. When COL
#' calibration maps different NCBI TaxIDs to the same accepted species, their
#' common taxonomy and \code{resolved_taxon_id} are identical and therefore do
#' not create an artificial conflict during LCA calculation.
#'
#' Assignments above species level receive a numeric label, for example
#' `Homo sp_1` and `Homo sp_2`. This keeps sequences or MOTUs with the same
#' higher-rank assignment distinguishable.
#'
#' Identity thresholds set only the deepest permitted rank. If identity permits
#' genus but the retained BLAST hits agree only at family, the result remains at
#' family.
#'
#' @examples
#' blast_hits <- data.frame(
#'   qseqid = c("query_1", "query_1"),
#'   sseqid = c("AB000001.1", "AB000002.1"),
#'   staxids = c(1, 2),
#'   pident = c(99.5, 99.5),
#'   qcovs = c(100, 100),
#'   evalue = c(1e-40, 1e-40),
#'   kingdom = c("Animalia", "Animalia"),
#'   phylum = c("Chordata", "Chordata"),
#'   class = c("Mammalia", "Mammalia"),
#'   order = c("Primates", "Primates"),
#'   family = c("Hominidae", "Hominidae"),
#'   genus = c("Homo", "Homo"),
#'   species = c("Homo sapiens", "Homo sapiens")
#' )
#' find_lca(blast_hits)
#' find_lca(
#'   blast_hits,
#'   identity_threshold = c(species = 100, genus = 98, family = 95)
#' )
#'
#' @importFrom dplyr across all_of coalesce filter group_by mutate select
#'   summarise sym syms
#' @export
find_lca <- function(
    df,
    seqs_ID = "qseqid",
    tax_levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
    taxid_column = "staxids",
    accession_column = "sseqid",
    identity_column = "pident",
    coverage_column = "qcovs",
    identity_threshold = NULL) {

  required <- c(
    seqs_ID, accession_column, taxid_column, identity_column,
    coverage_column, "evalue", tax_levels
  )
  missing_columns <- setdiff(required, names(df))
  if (length(missing_columns)) {
    stop(
      "Annotated BLAST data are missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  get_consensus <- function(x) {
    values <- unique(stats::na.omit(x))
    if (length(values) == 1L) values else NA_character_
  }

  best_hits <- df %>%
    group_by(!!sym(seqs_ID)) %>%
    filter(
      !!sym("evalue") == min(!!sym("evalue"), na.rm = TRUE)
    ) %>%
    ungroup()

  result <- best_hits %>%
    group_by(!!sym(seqs_ID)) %>%
    summarise(
      across(all_of(tax_levels), get_consensus),
      taxid = paste(unique(stats::na.omit(!!sym(taxid_column))), collapse = ";"),
      accession = paste(unique(stats::na.omit(!!sym(accession_column))), collapse = ";"),
      coverage = unique(!!sym(coverage_column))[1],
      identity = unique(!!sym(identity_column))[1],
      .groups = "drop"
    )

  if (!is.null(identity_threshold)) {
    threshold_ranks <- names(identity_threshold)
    identity_threshold <- as.numeric(identity_threshold)
    names(identity_threshold) <- threshold_ranks
    unknown_ranks <- setdiff(threshold_ranks, tax_levels)
    if (is.null(threshold_ranks) || length(unknown_ranks)) {
      stop(
        "identity_threshold must be named with ranks from tax_levels.",
        call. = FALSE
      )
    }

    for (rank in threshold_ranks) {
      rank_position <- match(rank, tax_levels)
      below_threshold <- is.na(result$identity) |
        result$identity < identity_threshold[[rank]]
      restricted_ranks <- tax_levels[rank_position:length(tax_levels)]
      result[below_threshold, restricted_ranks] <- NA_character_
    }
  }

  if ("resolved" %in% names(best_hits)) {
    query_values <- unique(best_hits[[seqs_ID]])
    resolution <- data.frame(query_values, stringsAsFactors = FALSE)
    names(resolution) <- seqs_ID
    resolution$resolved <- vapply(query_values, function(query) {
      rows <- best_hits[[seqs_ID]] == query
      all(best_hits$resolved[rows] %in% TRUE)
    }, logical(1))
    resolution$resolved_taxon_id <- if ("resolved_taxon_id" %in% names(best_hits)) {
      vapply(query_values, function(query) {
        rows <- best_hits[[seqs_ID]] == query
        ids <- unique(stats::na.omit(best_hits$resolved_taxon_id[rows]))
        if (length(ids)) paste(ids, collapse = ";") else NA_character_
      }, character(1))
    } else {
      NA_character_
    }
    result <- dplyr::left_join(result, resolution, by = seqs_ID)
  }

  result$rank <- apply(
    result[, rev(tax_levels), drop = FALSE],
    1,
    function(x) {
      known <- which(!is.na(x))
      if (length(known)) rev(tax_levels)[known[1]] else NA_character_
    }
  )
  result$scientific_name <- apply(
    result[, rev(tax_levels), drop = FALSE],
    1,
    function(x) {
      values <- x[!is.na(x)]
      if (length(values)) values[1] else NA_character_
    }
  )
  above_species <- !is.na(result$rank) &
    result$rank != tax_levels[length(tax_levels)] &
    !is.na(result$scientific_name)
  result$scientific_name[above_species] <- paste0(
    result$scientific_name[above_species],
    " sp_",
    seq_len(sum(above_species))
  )
  result
}
