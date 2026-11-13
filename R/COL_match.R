#' Build a Catalogue of Life taxonomy database
#'
#' Build an indexed SQLite database from either a Catalogue of Life
#' `NameUsage.tsv` path or an in-memory data frame. File input is read and
#' written in chunks, which is the recommended memory-friendly method for a
#' complete COL export. In-memory input can construct missing lineage columns
#' from parent relationships.
#'
#' The input must contain Catalogue of Life-style columns such as `col:ID`,
#' `col:parentID`, `col:scientificName`, `col:rank`, and `col:status`. Column
#' names may retain the `col:` prefix or use plain names. In ColDP `NameUsage`
#' data, `parentID` identifies the taxonomic parent for accepted taxa and the
#' accepted taxon for synonyms. Explicit `acceptedNameUsageID` or
#' `acceptedTaxonID` fields take precedence when present.
#'
#' @param col_data Path to a COL `NameUsage.tsv` file, or a data frame or
#'   tibble containing COL taxon data.
#' @param db_path Path to the SQLite database to create.
#' @param table Name of the SQLite table to write.
#' @param overwrite Logical. If `TRUE`, remove an existing database at
#'   `db_path`.
#' @param max_depth Maximum parent-chain iterations used when filling lineage
#'   columns for in-memory input.
#' @param chunk_size Number of rows per read or write chunk.
#'
#' @return Invisibly returns `db_path`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Recommended for a complete COL export
#' build_col_taxonomy_database("NameUsage.tsv", "col.sqlite", overwrite = TRUE)
#'
#' # An already loaded COL table is also supported
#' build_col_taxonomy_database(col_data, "col.sqlite", overwrite = TRUE)
#' }
build_col_taxonomy_database <- function(col_data,
                             db_path,
                             table = "col_taxa",
                             overwrite = FALSE,
                             max_depth = 40,
                             chunk_size = 100000L) {
  if (is.character(col_data) && length(col_data) != 1L) {
    stop(
      "Character input to col_data must be one NameUsage.tsv file path.",
      call. = FALSE
    )
  }
  is_tsv_path <- is.character(col_data) && length(col_data) == 1L
  if (!is_tsv_path) {
    return(build_col_taxonomy_database_in_memory(
      col_data = col_data,
      db_path = db_path,
      table = table,
      overwrite = overwrite,
      max_depth = max_depth,
      chunk_size = chunk_size
    ))
  }

  tsv_path <- path.expand(col_data)
  if (!file.exists(tsv_path)) {
    stop("COL TSV file does not exist: ", tsv_path, call. = FALSE)
  }

  col_require_namespace("readr")
  col_require_namespace("dplyr")
  col_require_namespace("DBI")
  col_require_namespace("RSQLite")

  header_con <- file(tsv_path, open = "r", encoding = "UTF-8")
  header_line <- tryCatch(
    readLines(header_con, n = 1L, warn = FALSE),
    finally = close(header_con)
  )

  header <- names(readr::read_tsv(
    I(header_line),
    n_max = 0,
    quote = "",
    show_col_types = FALSE
  ))
  clean_header <- col_clean_names(header)
  keep_clean <- col_required_source_columns(clean_header)
  keep_original <- header[clean_header %in% keep_clean]
  keep_col_types <- do.call(
    readr::cols_only,
    stats::setNames(rep(list(readr::col_character()), length(keep_original)), keep_original)
  )

  lineage_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  if (!all(lineage_ranks %in% keep_clean)) {
    stop(
      "Chunked TSV import requires COL lineage columns: ",
      paste(lineage_ranks, collapse = ", "),
      ". Use build_col_taxonomy_database() on an in-memory table if lineage must be computed from parentID.",
      call. = FALSE
    )
  }

  if (file.exists(db_path)) {
    if (!isTRUE(overwrite)) {
      stop("Database already exists. Use overwrite = TRUE to replace it.", call. = FALSE)
    }
    unlink(db_path)
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA synchronous = OFF")
  DBI::dbExecute(con, "PRAGMA temp_store = MEMORY")

  table_initialized <- FALSE
  callback <- function(x, pos) {
    taxa <- col_prepare_col_taxa(x, max_depth = max_depth, compute_lineage = FALSE)
    if (!table_initialized) {
      if (DBI::dbExistsTable(con, table)) {
        DBI::dbRemoveTable(con, table)
      }
      DBI::dbWriteTable(con, table, taxa[0, , drop = FALSE], overwrite = FALSE)
      table_initialized <<- TRUE
    }
    if (nrow(taxa) > 0L) {
      DBI::dbAppendTable(con, table, taxa)
    }
    invisible()
  }

  progress_options <- options(width = 50)
  on.exit(options(progress_options), add = TRUE)
  readr::read_tsv_chunked(
    file = tsv_path,
    callback = readr::SideEffectChunkCallback$new(callback),
    chunk_size = as.integer(chunk_size),
    col_types = keep_col_types,
    quote = "",
    progress = TRUE,
    show_col_types = FALSE
  )

  if (!table_initialized) {
    stop("No data rows were read from TSV.", call. = FALSE)
  }

  col_create_sqlite_indexes(con, table)
  DBI::dbExecute(con, sprintf("ANALYZE %s", DBI::dbQuoteIdentifier(con, table)))

  invisible(db_path)
}

build_col_taxonomy_database_in_memory <- function(col_data,
                                       db_path,
                                       table,
                                       overwrite,
                                       max_depth,
                                       chunk_size) {
  col_require_namespace("DBI")
  col_require_namespace("RSQLite")
  col_require_namespace("dplyr")
  col_require_namespace("tibble")

  if (file.exists(db_path)) {
    if (!isTRUE(overwrite)) {
      stop("Database already exists. Use overwrite = TRUE to replace it.", call. = FALSE)
    }
    unlink(db_path)
  }

  taxa <- col_prepare_col_taxa(col_data, max_depth = max_depth, compute_lineage = TRUE)

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA synchronous = OFF")
  DBI::dbExecute(con, "PRAGMA temp_store = MEMORY")

  if (DBI::dbExistsTable(con, table)) {
    DBI::dbRemoveTable(con, table)
  }
  DBI::dbWriteTable(con, table, taxa[0, , drop = FALSE], overwrite = FALSE)
  row_count <- nrow(taxa)
  if (row_count > 0L) {
    starts <- seq.int(1L, row_count, by = as.integer(chunk_size))
    for (start in starts) {
      end <- min(start + as.integer(chunk_size) - 1L, row_count)
      DBI::dbAppendTable(con, table, taxa[start:end, , drop = FALSE])
    }
  }

  col_create_sqlite_indexes(con, table)
  DBI::dbExecute(con, sprintf("ANALYZE %s", DBI::dbQuoteIdentifier(con, table)))

  invisible(db_path)
}

#' Match species names against a COL SQLite database
#'
#' Match a vector or table of species names to a SQLite database created by
#' `build_col_taxonomy_database()`. Results include the matched taxon record plus the full
#' standard lineage: `kingdom`, `phylum`, `class`, `order`, `family`, `genus`,
#' and `species`.
#'
#' Input names must not include taxonomic authorship. For example, use
#' `"Quercus robur"` rather than `"Quercus robur L."`. Authorship is stored
#' separately in the COL database and is not removed from queries by this
#' function. Fuzzy matching corrects limited spelling differences but does not
#' make authorship-containing names valid inputs.
#'
#' Exact matching is performed first using normalized scientific names. When
#' `fuzzy = TRUE`, unmatched names are compared against constrained candidate
#' sets using edit distance. The fuzzy normalization can treat common OCR or
#' typing confusions such as `i` and `l` as similar by using
#' `confusables = TRUE`.
#' For fuzzy matches, \code{match_distance} is the distance actually used for
#' acceptance: the smaller of the ordinary normalized-name distance and the
#' confusable-normalized distance.
#'
#' Constraints can be supplied in two ways. First, if `taxon_names` is a data frame and
#' contains columns named like taxonomy ranks, for example `kingdom` or
#' `family`, each row is matched only within those values. Second, fixed
#' constraints can be passed through arguments such as `kingdom = "Animalia"`.
#'
#' @param taxon_names A character vector of species names without authorship, or a data
#'   frame containing such a species-name column.
#' @param db_path Path to a SQLite database created by `build_col_taxonomy_database()`.
#' @param name_col Name of the species-name column when `taxon_names` is a data frame.
#' @param table SQLite table name used by `build_col_taxonomy_database()`.
#' @param ranks Taxonomic ranks eligible for matching. Defaults to `"species"`.
#' @param kingdom,phylum,class,order,family,genus Optional fixed constraints.
#'   `NULL` means no fixed constraint for that rank.
#' @param fuzzy Logical. If `TRUE`, run fuzzy matching for names not resolved by
#'   exact matching.
#' @param max_dist Maximum edit distance accepted for fuzzy matching.
#' @param method String distance method used by the optional `stringdist`
#'   package. If `stringdist` is not installed, base R `adist()` is used.
#' @param confusables Logical. If `TRUE`, the fuzzy key normalizes visually or
#'   typographically similar letters, currently `l` to `i`.
#' @param return_all Logical. If `TRUE`, return all exact matches. If `FALSE`,
#'   keep the first best match per input row, preferring full exact-name matches,
#'   accepted names, and lower fuzzy distance.
#' @param select Extra columns to return from SQLite. Core match columns and
#'   lineage columns are always returned. `NULL` returns all columns.
#' @param parallel Logical. If `TRUE`, use `future.apply` when enough unmatched
#'   names remain to offset the cost of starting parallel workers.
#' @param n_threads Number of parallel workers used for fuzzy matching. Defaults
#'   to 4.
#'
#' @return A tibble with original input columns prefixed by `input_` and matched
#'   COL columns. Matching results include `match_type`,
#'   `match_distance`, `input_name_status`, `matched_name_status`, `taxon_id`,
#'   `scientific_name`, `rank`, `kingdom`, `phylum`, `class`, `order`,
#'   `family`, `genus`, `species`, `accepted_name`, and
#'   `accepted_authorship`. When the input name is a synonym,
#'   `input_name_status` records `"synonym"`, while
#'   `matched_name_status` and the returned lineage describe its accepted
#'   taxon.
#' @export
#'
#' @examples
#' \dontrun{
#' col_match_species(c("Homo sapiens", "Panthera leo"), "col.sqlite")
#'
#' input <- tibble::tibble(
#'   species = c("Homo saplens", "Panthera ieo"),
#'   kingdom = c("Animalia", "Animalia")
#' )
#' col_match_species(input, "col.sqlite", fuzzy = TRUE, max_dist = 1)
#' }
col_match_species <- function(taxon_names,
                              db_path,
                              name_col = "species",
                              table = "col_taxa",
                              ranks = col_species_match_ranks(),
                              kingdom = NULL,
                              phylum = NULL,
                              class = NULL,
                              order = NULL,
                              family = NULL,
                              genus = NULL,
                              fuzzy = FALSE,
                              max_dist = 1L,
                              method = "lv",
                              confusables = TRUE,
                              return_all = FALSE,
                              select = NULL,
                              parallel = TRUE,
                              n_threads = 4L) {
  col_require_namespace("DBI")
  col_require_namespace("RSQLite")
  col_require_namespace("dplyr")
  col_require_namespace("tibble")

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  fields <- DBI::dbListFields(con, table)
  constraints <- col_collect_constraints(
    kingdom = kingdom,
    phylum = phylum,
    class = class,
    order = order,
    family = family,
    genus = genus
  )

  name_col <- if (is.character(taxon_names)) "species" else col_clean_names(name_col)
  input <- col_prepare_input(taxon_names, name_col)
  rank_cols <- intersect(c("kingdom", "phylum", "class", "order", "family", "genus"), names(input))
  for (col in names(constraints)) {
    input[[col]] <- constraints[[col]]
  }
  rank_cols <- unique(c(rank_cols, names(constraints)))

  input <- dplyr::mutate(
    input,
    input_row = dplyr::row_number(),
    query_name = as.character(.data[[name_col]]),
    query_norm = col_normalize_name(.data$query_name),
    query_fuzzy = col_normalize_name(.data$query_name, confusables = confusables)
  )
  dedupe_cols <- unique(c("query_norm", "query_fuzzy", rank_cols))
  lookup_map <- input %>%
    dplyr::select(dplyr::all_of(c("input_row", dedupe_cols))) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(dedupe_cols)), .keep_all = TRUE) %>%
    dplyr::mutate(lookup_row = dplyr::row_number())
  input <- dplyr::left_join(
    input,
    dplyr::select(lookup_map, dplyr::all_of(c(dedupe_cols, "lookup_row"))),
    by = dedupe_cols
  )
  lookup_input <- lookup_map %>%
    dplyr::select(-dplyr::all_of("input_row")) %>%
    dplyr::rename(input_row = "lookup_row")

  selected <- col_select_fields(fields, select)
  ranks <- col_rank_query_values(ranks)

  exact <- col_exact_matches(
    con = con,
    table = table,
    input = lookup_input,
    selected = selected,
    ranks = ranks,
    rank_cols = rank_cols
  )

  if (isTRUE(fuzzy)) {
    matched_rows <- if ("input_row" %in% names(exact)) unique(exact$input_row) else integer()
    todo <- dplyr::filter(
      lookup_input,
      !.data$input_row %in% matched_rows,
      !is.na(.data$query_fuzzy),
      nzchar(.data$query_fuzzy)
    )
    fuzzy_hits <- col_fuzzy_matches(
      con = con,
      db_path = db_path,
      table = table,
      input = todo,
      selected = selected,
      ranks = ranks,
      rank_cols = rank_cols,
      max_dist = as.integer(max_dist),
      method = method,
      parallel = parallel,
      n_threads = n_threads
    )
    out <- dplyr::bind_rows(exact, fuzzy_hits)
  } else {
    out <- exact
  }

  input_keep <- setdiff(names(input), c("query_norm", "query_fuzzy", "query_name"))
  input_part <- dplyr::select(input, dplyr::all_of(input_keep))
  names(input_part)[!names(input_part) %in% c("input_row", "lookup_row")] <- paste0(
    "input_",
    names(input_part)[!names(input_part) %in% c("input_row", "lookup_row")]
  )

  if (nrow(out) == 0L) {
    return(tibble::as_tibble(input_part[0, !names(input_part) %in% c("input_row", "lookup_row"), drop = FALSE]))
  }

  out <- col_add_accepted_names(con, table, out)
  out$full_exact_sort <- out$match_type == "exact" & out$scientific_name_norm == lookup_input$query_norm[match(out$input_row, lookup_input$input_row)]
  out$accepted_sort <- if ("status" %in% names(out)) out$status == "accepted" else FALSE
  out <- dplyr::arrange(
    out,
    .data$input_row,
    dplyr::desc(.data$full_exact_sort),
    .data$match_distance,
    dplyr::desc(.data$accepted_sort)
  )
  if (!isTRUE(return_all)) {
    out <- dplyr::slice(dplyr::group_by(out, .data$input_row), 1L)
    out <- dplyr::ungroup(out)
  }
  out <- dplyr::select(out, -dplyr::all_of(c("accepted_sort", "full_exact_sort")))
  out <- dplyr::rename(out, lookup_row = "input_row")
  out <- col_drop_internal_columns(out)

  dplyr::inner_join(input_part, out, by = "lookup_row") %>%
    dplyr::arrange(.data$input_row) %>%
    dplyr::select(-dplyr::all_of(c("input_row", "lookup_row"))) %>%
    tibble::as_tibble()
}

#' Cross-reference taxon names between Catalogue of Life and NCBI
#'
#' Standardize uncorrected species names, remove authorship, expand every COL
#' match to its accepted name and all synonyms, and then look up every expanded
#' name in a local NCBI taxonomy SQLite database. Inputs without a COL match are
#' looked up directly in NCBI using the same standardized name used for COL.
#' Infraspecific markers
#' such as `subsp.`, `var.`, `f.`, `subvar.`, and `subf.` are retained and
#' normalized.
#'
#' @param taxon_names Character vector containing uncorrected scientific names.
#' @param col_db Path to the SQLite database created by `build_col_taxonomy_database()`.
#' @param ncbi_db Path to the NCBI taxonomy SQLite database created by
#'   `build_ncbi_taxonomy_database(taxonomy = TRUE)`.
#' @param col_table Name of the COL table. Default is `"col_taxa"`.
#' @param ranks,kingdom,phylum,class,order,family,genus,fuzzy,max_dist,method,
#'   confusables,return_all,select,parallel,n_threads Matching options passed to
#'   `col_match_species()`.
#'
#' @return A tibble containing the original and standardized input names, the
#'   COL name initially matched, its accepted taxon and accepted name, every COL
#'   accepted/synonym name in that concept, and corresponding NCBI name classes
#'   and taxids. For a direct NCBI fallback match, the COL columns are missing.
#'   Inputs absent from both databases and COL names without an NCBI match are
#'   retained with missing values.
#' @export
#'
#' @examples
#' \dontrun{
#' col_ncbi_crosswalk(
#'   c("Quercus robur L.", "Acacia nilotica (L.) Delile subsp. indica"),
#'   col_db = "col.sqlite",
#'   ncbi_db = "ncbi_taxonomy.sqlite",
#'   fuzzy = TRUE,
#'   max_dist = 1
#' )
#' }
col_ncbi_crosswalk <- function(taxon_names,
                               col_db,
                               ncbi_db,
                               col_table = "col_taxa",
                               ranks = col_species_match_ranks(),
                               kingdom = NULL,
                               phylum = NULL,
                               class = NULL,
                               order = NULL,
                               family = NULL,
                               genus = NULL,
                               fuzzy = FALSE,
                               max_dist = 1L,
                               method = "lv",
                               confusables = TRUE,
                               return_all = FALSE,
                               select = NULL,
                               parallel = TRUE,
                               n_threads = 4L) {
  input_name <- trimws(gsub("[[:space:]]+", " ", as.character(taxon_names)))
  standardized_name <- vapply(input_name, function(x) {
    words <- regmatches(x, gregexpr("[[:alpha:]][[:alpha:]-]*", x))[[1]]
    if (!length(words)) return(NA_character_)
    if (length(words) == 1L) return(paste0(toupper(substr(words[1], 1, 1)), tolower(substr(words[1], 2, nchar(words[1])))))

    genus <- paste0(toupper(substr(words[1], 1, 1)), tolower(substr(words[1], 2, nchar(words[1]))))
    species <- tolower(words[2])
    lower_words <- tolower(words)
    rank_map <- c(
      subsp = "subsp.", ssp = "subsp.", subspecies = "subsp.",
      var = "var.", varietas = "var.", variety = "var.",
      f = "f.", fo = "f.", forma = "f.", form = "f.",
      subvar = "subvar.", subvarietas = "subvar.", subvariety = "subvar.",
      subf = "subf.", subforma = "subf.", subform = "subf."
    )
    rank_position <- which(lower_words %in% names(rank_map) & seq_along(words) > 2L)
    if (length(rank_position) && rank_position[1] < length(words)) {
      i <- rank_position[1]
      return(paste(genus, species, rank_map[lower_words[i]], tolower(words[i + 1L])))
    }

    author_connectors <- c("ex", "in", "non", "sensu", "auct", "apud",
                           "de", "del", "della", "van", "von", "der", "den", "du", "la", "le")
    third_is_epithet <- length(words) >= 3L &&
      grepl("^[[:lower:]][[:alpha:]-]*$", words[3]) &&
      !lower_words[3] %in% author_connectors
    if (third_is_epithet) paste(genus, species, tolower(words[3])) else paste(genus, species)
  }, character(1), USE.NAMES = FALSE)

  input <- data.frame(
    input_id = seq_along(input_name),
    input_name = input_name,
    standardized_name = standardized_name,
    query_norm = tolower(gsub("[^[:alnum:]]+", " ", standardized_name)),
    stringsAsFactors = FALSE
  )
  input$query_norm <- trimws(gsub("[[:space:]]+", " ", input$query_norm))

  col_con <- DBI::dbConnect(RSQLite::SQLite(), col_db)
  on.exit(DBI::dbDisconnect(col_con), add = TRUE)
  ncbi_con <- DBI::dbConnect(RSQLite::SQLite(), ncbi_db)
  on.exit(DBI::dbDisconnect(ncbi_con), add = TRUE)

  ncbi_table <- if ("name_aliases" %in% DBI::dbListTables(ncbi_con)) "name_aliases" else "names"
  query_ncbi_names <- function(search_names) {
    search_names <- unique(search_names[!is.na(search_names) & nzchar(search_names)])
    if (!length(search_names)) return(data.frame())

    hits <- lapply(split(search_names, ceiling(seq_along(search_names) / 5000L)), function(chunk) {
      if (ncbi_table == "name_aliases") {
        sql <- sprintf(
          "SELECT tax_id, name, name_class FROM name_aliases WHERE name IN (%s)",
          paste(DBI::dbQuoteString(ncbi_con, chunk), collapse = ", ")
        )
      } else {
        sql <- sprintf(
          "SELECT tax_id, name, 'scientific name' AS name_class FROM names WHERE name IN (%s)",
          paste(DBI::dbQuoteString(ncbi_con, chunk), collapse = ", ")
        )
      }
      DBI::dbGetQuery(ncbi_con, sql)
    })
    hits <- do.call(rbind, hits)
    if (!nrow(hits)) return(data.frame())

    hits$tax_id <- .resolve_ncbi_taxids(ncbi_con, hits$tax_id)
    current_names <- DBI::dbGetQuery(
      ncbi_con,
      sprintf(
        "SELECT tax_id, name FROM names WHERE tax_id IN (%s)",
        paste(DBI::dbQuoteLiteral(ncbi_con, unique(hits$tax_id)), collapse = ",")
      )
    )
    names(current_names) <- c("ncbi_taxid", "ncbi_scientific_name")
    names(hits) <- c("ncbi_taxid", "ncbi_matched_name", "ncbi_name_class")
    merge(hits, current_names, by = "ncbi_taxid", all.x = TRUE)
  }

  make_ncbi_fallback <- function(unmatched_input) {
    # Preserve the supplied name in input_name, but query NCBI with the same
    # authorship-stripped, normalized name used for the COL lookup.
    hits <- query_ncbi_names(unmatched_input$standardized_name)
    if (nrow(hits)) {
      fallback <- merge(
        unmatched_input[, c("input_id", "input_name", "standardized_name"), drop = FALSE],
        hits,
        by.x = "standardized_name", by.y = "ncbi_matched_name",
        all.x = TRUE, sort = FALSE
      )
    } else {
      fallback <- unmatched_input[, c("input_id", "input_name", "standardized_name"), drop = FALSE]
      fallback$ncbi_taxid <- NA_integer_
      fallback$ncbi_name_class <- NA_character_
      fallback$ncbi_scientific_name <- NA_character_
    }
    fallback$accepted_taxon_id <- NA_character_
    fallback$col_input_match_name <- NA_character_
    fallback$col_input_match_status <- NA_character_
    fallback$col_accepted_name <- NA_character_
    fallback$col_taxon_id <- NA_character_
    fallback$col_name <- NA_character_
    fallback$col_name_status <- NA_character_
    fallback$col_rank <- NA_character_
    fallback$col_authorship <- NA_character_
    fallback
  }

  matches <- col_match_species(
    taxon_names = data.frame(
      crosswalk_id = input$input_id,
      species = input$standardized_name,
      stringsAsFactors = FALSE
    ),
    db_path = col_db,
    name_col = "species",
    table = col_table,
    ranks = ranks,
    kingdom = kingdom,
    phylum = phylum,
    class = class,
    order = order,
    family = family,
    genus = genus,
    fuzzy = fuzzy,
    max_dist = max_dist,
    method = method,
    confusables = confusables,
    return_all = return_all,
    select = select,
    parallel = parallel,
    n_threads = n_threads
  )

  initial <- data.frame()
  if (nrow(matches)) {
    initial <- data.frame(
      input_id = matches$input_crosswalk_id,
      col_input_match_name = matches$scientific_name,
      col_input_match_status = matches$input_name_status,
      accepted_taxon_id = matches$accepted_taxon_id,
      match_type = matches$match_type,
      match_distance = matches$match_distance,
      stringsAsFactors = FALSE
    )
  }

  if (!nrow(initial)) {
    result <- make_ncbi_fallback(input)
    result <- result[order(result$input_id), , drop = FALSE]
    rownames(result) <- NULL
    return(tibble::as_tibble(result))
  }

  accepted_ids <- unique(initial$accepted_taxon_id[!is.na(initial$accepted_taxon_id)])
  col_fields <- DBI::dbListFields(col_con, col_table)
  wanted <- intersect(
    c("taxon_id", "parent_id", "accepted_name_usage_id", "scientific_name",
      "status", "rank", "authorship"),
    col_fields
  )
  expanded <- lapply(split(accepted_ids, ceiling(seq_along(accepted_ids) / 5000L)), function(chunk) {
    ids <- paste(DBI::dbQuoteString(col_con, chunk), collapse = ", ")
    sql <- sprintf(
      paste0(
        "SELECT %s FROM %s WHERE taxon_id IN (%s) ",
        "OR accepted_name_usage_id IN (%s) ",
        "OR (lower(status) IN ('synonym','ambiguous synonym','misapplied') AND parent_id IN (%s))"
      ),
      paste(DBI::dbQuoteIdentifier(col_con, wanted), collapse = ", "),
      DBI::dbQuoteIdentifier(col_con, col_table), ids, ids, ids
    )
    DBI::dbGetQuery(col_con, sql)
  })
  expanded <- do.call(rbind, expanded)
  synonym_status <- c("synonym", "ambiguous synonym", "misapplied")
  for (column in c("parent_id", "accepted_name_usage_id", "status", "rank", "authorship")) {
    if (!column %in% names(expanded)) expanded[[column]] <- NA_character_
  }
  expanded_status <- tolower(expanded$status)
  expanded$accepted_taxon_id <- as.character(expanded$taxon_id)
  expanded_is_synonym <- expanded_status %in% synonym_status
  has_parent <- !is.na(expanded$parent_id) & nzchar(as.character(expanded$parent_id))
  has_explicit <- !is.na(expanded$accepted_name_usage_id) & nzchar(as.character(expanded$accepted_name_usage_id))
  expanded$accepted_taxon_id[expanded_is_synonym & has_parent] <- as.character(expanded$parent_id[expanded_is_synonym & has_parent])
  expanded$accepted_taxon_id[expanded_is_synonym & has_explicit] <- as.character(expanded$accepted_name_usage_id[expanded_is_synonym & has_explicit])

  accepted <- expanded[as.character(expanded$taxon_id) == expanded$accepted_taxon_id,
                       c("taxon_id", "scientific_name"), drop = FALSE]
  accepted <- accepted[!duplicated(accepted$taxon_id), , drop = FALSE]
  names(accepted) <- c("accepted_taxon_id", "col_accepted_name")
  expanded <- merge(expanded, accepted, by = "accepted_taxon_id", all.x = TRUE)
  names(expanded)[names(expanded) == "taxon_id"] <- "col_taxon_id"
  names(expanded)[names(expanded) == "scientific_name"] <- "col_name"
  names(expanded)[names(expanded) == "status"] <- "col_name_status"
  names(expanded)[names(expanded) == "rank"] <- "col_rank"
  names(expanded)[names(expanded) == "authorship"] <- "col_authorship"

  crosswalk <- merge(
    initial[, c("input_id", "col_input_match_name", "col_input_match_status",
                "accepted_taxon_id", "match_type", "match_distance"), drop = FALSE],
    expanded[, c("accepted_taxon_id", "col_accepted_name", "col_taxon_id", "col_name",
                 "col_name_status", "col_rank", "col_authorship"), drop = FALSE],
    by = "accepted_taxon_id"
  )

  ncbi_names <- unique(crosswalk$col_name)
  ncbi_hits <- query_ncbi_names(ncbi_names)

  if (nrow(ncbi_hits)) {
    crosswalk <- merge(crosswalk, ncbi_hits, by.x = "col_name", by.y = "ncbi_matched_name", all.x = TRUE)
  } else {
    crosswalk$ncbi_taxid <- NA_integer_
    crosswalk$ncbi_name_class <- NA_character_
    crosswalk$ncbi_scientific_name <- NA_character_
  }

  result <- merge(
    input[, c("input_id", "input_name", "standardized_name")],
    crosswalk[, setdiff(names(crosswalk), c("input_name", "standardized_name")), drop = FALSE],
    by = "input_id", all.x = TRUE, sort = FALSE
  )

  unmatched_input <- input[!input$input_id %in% initial$input_id, , drop = FALSE]
  if (nrow(unmatched_input)) {
    result <- result[result$input_id %in% initial$input_id, , drop = FALSE]
    fallback <- make_ncbi_fallback(unmatched_input)
    all_columns <- union(names(result), names(fallback))
    for (column in setdiff(all_columns, names(result))) result[[column]] <- NA
    for (column in setdiff(all_columns, names(fallback))) fallback[[column]] <- NA
    result <- rbind(result[, all_columns, drop = FALSE], fallback[, all_columns, drop = FALSE])
  }
  result <- result[order(result$input_id, result$accepted_taxon_id, result$col_name), , drop = FALSE]
  rownames(result) <- NULL
  tibble::as_tibble(result)
}

col_require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required. Install it first.", call. = FALSE)
  }
}

col_required_source_columns <- function(names_clean) {
  required <- c("id", "parentid", "scientificname", "rank")
  optional <- c(
    "alternativeid", "namealternativeid", "sourceid", "basionymid",
    "acceptednameusageid", "acceptedtaxonid", "status",
    "authorship", "notho", "originalspelling", "uninomial", "genericname",
    "infragenericepithet", "specificepithet", "infraspecificepithet",
    "cultivarepithet", "combinationauthorship",
    "kingdom", "phylum", "class", "order", "family", "genus", "species"
  )
  unique(c(required, optional[optional %in% names_clean]))
}

col_prepare_col_taxa <- function(col_data, max_depth, compute_lineage) {
  col_require_namespace("dplyr")
  col_require_namespace("tibble")

  taxa <- tibble::as_tibble(col_data)
  names(taxa) <- col_clean_names(names(taxa))

  required <- c("id", "parentid", "scientificname", "rank")
  missing <- setdiff(required, names(taxa))
  if (length(missing) > 0L) {
    stop("Missing required COL columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  keep <- col_required_source_columns(names(taxa))
  taxa <- dplyr::select(taxa, dplyr::all_of(keep))

  rename_map <- c(
    id = "taxon_id",
    parentid = "parent_id",
    scientificname = "scientific_name",
    genericname = "generic_name",
    specificepithet = "specific_epithet",
    infraspecificepithet = "infraspecific_epithet"
  )
  for (old in intersect(names(rename_map), names(taxa))) {
    names(taxa)[names(taxa) == old] <- rename_map[[old]]
  }

  accepted_name_usage_id <- rep(NA_character_, nrow(taxa))
  if ("acceptednameusageid" %in% names(taxa)) {
    accepted_name_usage_id <- as.character(taxa$acceptednameusageid)
  }
  if ("acceptedtaxonid" %in% names(taxa)) {
    accepted_name_usage_id <- dplyr::coalesce(
      accepted_name_usage_id,
      as.character(taxa$acceptedtaxonid)
    )
  }
  taxa$accepted_name_usage_id <- accepted_name_usage_id
  taxa <- dplyr::select(
    taxa,
    -dplyr::any_of(c("acceptednameusageid", "acceptedtaxonid"))
  )
  if (!"status" %in% names(taxa)) {
    taxa$status <- NA_character_
  }

  taxa <- dplyr::mutate(
    taxa,
    dplyr::across(
      dplyr::everything(),
      ~ dplyr::na_if(trimws(as.character(.x)), "")
    ),
    status = tolower(.data$status),
    rank = col_normalize_rank(.data$rank),
    scientific_name_norm = col_normalize_name(.data$scientific_name),
    scientific_name_fuzzy = col_normalize_name(.data$scientific_name, confusables = TRUE)
  )

  lineage_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  missing_lineage <- setdiff(lineage_ranks, names(taxa))
  if (!isTRUE(compute_lineage)) {
    if (length(missing_lineage) > 0L) {
      stop(
        "Chunked TSV import requires existing lineage columns: ",
        paste(lineage_ranks, collapse = ", "),
        call. = FALSE
      )
    }
  } else {
    lineage <- tibble::as_tibble(col_make_lineage(
      taxon_id = taxa$taxon_id,
      parent_id = taxa$parent_id,
      rank = taxa$rank,
      scientific_name = taxa$scientific_name,
      ranks = lineage_ranks,
      max_depth = max_depth
    ))
    for (rank_name in lineage_ranks) {
      if (rank_name %in% names(taxa)) {
        taxa[[rank_name]] <- dplyr::coalesce(taxa[[rank_name]], lineage[[rank_name]])
      } else {
        taxa[[rank_name]] <- lineage[[rank_name]]
      }
    }
  }
  for (rank_name in lineage_ranks) {
    taxa[[paste0(rank_name, "_norm")]] <- col_normalize_name(taxa[[rank_name]])
  }
  is_infraspecific <- taxa$rank %in% col_infraspecific_ranks()
  taxa$species[is_infraspecific] <- taxa$scientific_name[is_infraspecific]
  taxa$species_norm <- col_normalize_name(taxa$species)

  taxa
}

col_create_sqlite_indexes <- function(con, table) {
  lineage_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  index_cols <- intersect(
    c(
      "scientific_name_norm", "scientific_name_fuzzy", "rank",
      lineage_ranks, paste0(lineage_ranks, "_norm"), "taxon_id", "parent_id",
      "accepted_name_usage_id"
    ),
    DBI::dbListFields(con, table)
  )
  for (col in index_cols) {
    sql <- sprintf(
      "CREATE INDEX IF NOT EXISTS %s ON %s (%s)",
      DBI::dbQuoteIdentifier(con, paste0("idx_", table, "_", col)),
      DBI::dbQuoteIdentifier(con, table),
      DBI::dbQuoteIdentifier(con, col)
    )
    DBI::dbExecute(con, sql)
  }
  invisible(index_cols)
}

col_clean_names <- function(x) {
  x <- sub("^col:", "", x, ignore.case = TRUE)
  x <- gsub("[^A-Za-z0-9_]+", "", x)
  tolower(x)
}

col_species_match_ranks <- function() {
  c("species", col_infraspecific_ranks())
}

col_infraspecific_ranks <- function() {
  c(
    "subspecies", "variety", "form",
    "subvariety", "subform",
    "cultivar", "strain"
  )
}

col_normalize_rank <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[._-]+", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  x <- trimws(x)
  dplyr::case_when(
    x %in% c("sp", "species") ~ "species",
    x %in% c("ssp", "subsp", "subspecies") ~ "subspecies",
    x %in% c("var", "varietas", "variety") ~ "variety",
    x %in% c("f", "fo", "forma", "form") ~ "form",
    x %in% c("subvar", "subvarietas", "subvariety") ~ "subvariety",
    x %in% c("subf", "subform", "subforma") ~ "subform",
    TRUE ~ x
  )
}

col_rank_query_values <- function(ranks) {
  ranks <- unique(c(as.character(ranks), col_normalize_rank(ranks)))
  rank_aliases <- list(
    species = c("sp", "species"),
    subspecies = c("ssp", "subsp", "subspecies"),
    variety = c("var", "varietas", "variety"),
    form = c("f", "fo", "forma", "form"),
    subvariety = c("subvar", "subvarietas", "subvariety"),
    subform = c("subf", "subforma", "subform")
  )
  canonical <- col_normalize_rank(ranks)
  aliases <- unlist(rank_aliases[intersect(canonical, names(rank_aliases))], use.names = FALSE)
  unique(c(ranks, aliases))
}

col_normalize_name <- function(x, confusables = FALSE) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]]+", " ", x)
  x <- gsub("[^a-z0-9 ._-]+", "", x)
  x <- gsub("[._-]+", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  if (isTRUE(confusables)) {
    x <- chartr("l", "i", x)
  }
  x
}

col_make_lineage <- function(taxon_id, parent_id, rank, scientific_name, ranks, max_depth) {
  row_count <- length(taxon_id)
  lineage <- matrix(NA_character_, nrow = row_count, ncol = length(ranks))
  colnames(lineage) <- ranks

  rank_index <- match(rank, ranks)
  has_rank <- !is.na(rank_index) & !is.na(scientific_name) & nzchar(scientific_name)
  lineage[cbind(which(has_rank), rank_index[has_rank])] <- scientific_name[has_rank]

  parent_index <- match(parent_id, taxon_id)
  parent_index[is.na(parent_index)] <- NA_integer_

  for (depth in seq_len(max_depth)) {
    changed <- FALSE
    has_parent <- !is.na(parent_index)
    for (j in seq_along(ranks)) {
      missing <- is.na(lineage[, j]) & has_parent
      if (!any(missing)) {
        next
      }
      values <- lineage[parent_index[missing], j]
      fill <- !is.na(values)
      if (any(fill)) {
        idx <- which(missing)[fill]
        lineage[idx, j] <- values[fill]
        changed <- TRUE
      }
    }
    if (!changed) {
      break
    }
  }

  lineage
}

col_prepare_input <- function(x, name_col) {
  if (is.character(x)) {
    tibble::tibble(species = x)
  } else {
    input <- tibble::as_tibble(x)
    names(input) <- col_clean_names(names(input))
    if (!name_col %in% names(input)) {
      stop("name_col was not found in x: ", name_col, call. = FALSE)
    }
    input
  }
}

col_collect_constraints <- function(...) {
  args <- list(...)
  args <- args[!vapply(args, is.null, logical(1))]
  if (length(args) == 0L) {
    return(list())
  }
  lens <- vapply(args, length, integer(1))
  if (any(lens != 1L)) {
    stop("Fixed taxonomy constraints must each be length 1.", call. = FALSE)
  }
  args
}

col_select_fields <- function(fields, select) {
  core <- c(
    "taxon_id", "parent_id", "scientific_name", "scientific_name_norm",
    "scientific_name_fuzzy", "status", "rank", "authorship",
    "accepted_name_usage_id",
    "kingdom", "phylum", "class", "order", "family", "genus", "species"
  )
  core <- core[core %in% fields]
  if (is.null(select)) {
    return(fields)
  }
  missing <- setdiff(select, fields)
  if (length(missing) > 0L) {
    stop("Selected columns are absent from SQLite table: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  unique(c(core, select))
}

col_exact_matches <- function(con, table, input, selected, ranks, rank_cols) {
  groups <- col_group_input(input, rank_cols)
  out <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    g <- groups[[i]]
    names_to_match <- unique(g$query_norm[!is.na(g$query_norm) & nzchar(g$query_norm)])
    if (length(names_to_match) == 0L) {
      next
    }
    chunks <- split(names_to_match, ceiling(seq_along(names_to_match) / 5000L))
    chunk_hits <- vector("list", length(chunks))
    for (chunk_id in seq_along(chunks)) {
      where <- c(
        col_sql_in(con, "scientific_name_norm", chunks[[chunk_id]]),
        col_sql_in(con, "rank", ranks),
        col_group_where(con, g[1, , drop = FALSE], rank_cols)
      )
      sql <- col_sql_select(con, table, selected, where)
      chunk_hits[[chunk_id]] <- tibble::as_tibble(DBI::dbGetQuery(con, sql))
    }
    hits <- dplyr::bind_rows(chunk_hits)
    if (nrow(hits) == 0L) {
      next
    }
    hits$match_key <- hits$scientific_name_norm
    query_keys <- dplyr::select(g, input_row, query_norm)
    hits <- dplyr::inner_join(hits, query_keys, by = c("match_key" = "query_norm"))
    hits$match_type <- "exact"
    hits$match_distance <- 0L
    hits$match_key <- NULL
    out[[i]] <- hits
  }
  dplyr::bind_rows(out)
}

col_fuzzy_matches <- function(con,
                              db_path,
                              table,
                              input,
                              selected,
                              ranks,
                              rank_cols,
                              max_dist,
                              method,
                              parallel,
                              n_threads) {
  if (nrow(input) == 0L) {
    return(tibble::tibble())
  }

  rows <- lapply(seq_len(nrow(input)), function(i) input[i, , drop = FALSE])
  use_parallel <- isTRUE(parallel) &&
    length(rows) >= max(20L, as.integer(n_threads) * 2L)
  if (use_parallel) {
    groups <- split(seq_along(rows),
                    rep_len(seq_len(min(as.integer(n_threads), length(rows))),
                            length(rows)))
    out <- col_future_lapply(groups, function(indices) {
      worker_con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
      on.exit(DBI::dbDisconnect(worker_con), add = TRUE)
      lapply(rows[indices], function(row_data) {
        col_fuzzy_match_one(
          con = worker_con,
          table = table,
          row_data = row_data,
          selected = selected,
          ranks = ranks,
          rank_cols = rank_cols,
          max_dist = max_dist,
          method = method
        )
      })
    }, n_threads = n_threads)
    out <- unlist(out, recursive = FALSE)
  } else {
    out <- lapply(rows, function(row_data) {
      col_fuzzy_match_one(
        con = con,
        table = table,
        row_data = row_data,
        selected = selected,
        ranks = ranks,
        rank_cols = rank_cols,
        max_dist = max_dist,
        method = method
      )
    })
  }
  dplyr::bind_rows(out)
}

col_fuzzy_match_one <- function(con,
                                table,
                                row_data,
                                selected,
                                ranks,
                                rank_cols,
                                max_dist,
                                method) {
  min_len <- nchar(row_data$query_fuzzy[[1]]) - max_dist
  max_len <- nchar(row_data$query_fuzzy[[1]]) + max_dist
  initial <- substr(row_data$query_fuzzy[[1]], 1L, 1L)
  fuzzy_column <- DBI::dbQuoteIdentifier(con, "scientific_name_fuzzy")
  prefix_start <- DBI::dbQuoteString(con, initial)
  prefix_end <- DBI::dbQuoteString(con, paste0(initial, "\U0010ffff"))
  where <- c(
    col_sql_in(con, "rank", ranks),
    sprintf("%s >= %s AND %s < %s",
            fuzzy_column, prefix_start, fuzzy_column, prefix_end),
    sprintf("length(scientific_name_fuzzy) BETWEEN %d AND %d", max(1L, min_len), max_len),
    col_group_where(con, row_data, rank_cols)
  )
  candidate_fields <- intersect(
    c("taxon_id", "scientific_name_norm", "scientific_name_fuzzy"),
    selected
  )
  sql <- col_sql_select(con, table, candidate_fields, where)
  candidates <- tibble::as_tibble(DBI::dbGetQuery(con, sql))
  if (nrow(candidates) == 0L) {
    return(NULL)
  }

  fuzzy_distances <- col_string_distance(row_data$query_fuzzy[[1]], candidates$scientific_name_fuzzy, method)
  raw_distances <- col_string_distance(row_data$query_norm[[1]], candidates$scientific_name_norm, method)
  keep_distance <- pmin(fuzzy_distances, raw_distances, na.rm = TRUE)
  keep <- which(!is.na(keep_distance) & keep_distance <= max_dist)
  if (length(keep) == 0L) {
    return(NULL)
  }
  matched <- candidates[keep, , drop = FALSE]
  hit_sql <- col_sql_select(
    con,
    table,
    selected,
    col_sql_in(con, "taxon_id", matched$taxon_id)
  )
  hit <- tibble::as_tibble(DBI::dbGetQuery(con, hit_sql))
  distance_index <- match(hit$taxon_id, matched$taxon_id)
  hit$input_row <- row_data$input_row[[1]]
  hit$match_type <- "fuzzy"
  hit$match_distance <- as.integer(keep_distance[keep][distance_index])
  hit
}

col_add_accepted_names <- function(con, table, hits) {
  if (nrow(hits) == 0L || !"taxon_id" %in% names(hits)) {
    return(hits)
  }

  status <- if ("status" %in% names(hits)) {
    tolower(trimws(as.character(hits$status)))
  } else {
    rep(NA_character_, nrow(hits))
  }
  hits$input_name_status <- status
  accepted_name_usage_id <- if ("accepted_name_usage_id" %in% names(hits)) {
    as.character(hits$accepted_name_usage_id)
  } else {
    rep(NA_character_, nrow(hits))
  }
  synonym <- !is.na(status) & status %in% c(
    "synonym", "ambiguous synonym", "misapplied"
  )
  explicit_accepted_id <- !is.na(accepted_name_usage_id) &
    nzchar(accepted_name_usage_id)
  coldp_accepted_id <- !is.na(hits$parent_id) & nzchar(hits$parent_id)
  unresolved_synonym <- synonym & (
    !explicit_accepted_id & !coldp_accepted_id
  )
  if (any(unresolved_synonym)) {
    warning(
      "Some synonym records contain neither an explicit accepted-name ID nor ",
      "a ColDP parentID, so their accepted taxa cannot be resolved.",
      call. = FALSE
    )
  }
  synonym_target_id <- dplyr::coalesce(
    ifelse(explicit_accepted_id, accepted_name_usage_id, NA_character_),
    ifelse(coldp_accepted_id, as.character(hits$parent_id), NA_character_)
  )
  target_id <- ifelse(
    synonym,
    synonym_target_id,
    hits$taxon_id
  )
  target_id <- as.character(target_id)
  ids <- unique(target_id[!is.na(target_id) & nzchar(target_id)])
  if (length(ids) == 0L) {
    hits$accepted_taxon_id <- NA_character_
    hits$accepted_name <- NA_character_
    hits$accepted_authorship <- NA_character_
    hits$accepted_status <- NA_character_
    return(hits)
  }

  fields <- DBI::dbListFields(con, table)
  accepted_cols <- intersect(
    c(
      "taxon_id", "scientific_name", "authorship", "status", "rank",
      "kingdom", "phylum", "class", "order", "family", "genus", "species"
    ),
    fields
  )
  chunks <- split(ids, ceiling(seq_along(ids) / 5000L))
  accepted <- lapply(chunks, function(chunk) {
    sql <- col_sql_select(
      con,
      table,
      selected = accepted_cols,
      where = col_sql_in(con, "taxon_id", chunk)
    )
    tibble::as_tibble(DBI::dbGetQuery(con, sql))
  }) %>%
    dplyr::bind_rows()

  if (nrow(accepted) == 0L) {
    warning(
      "Accepted taxon IDs referenced by matched names were not found in the ",
      "COL SQLite database. Rebuild it from the complete NameUsage.tsv file.",
      call. = FALSE
    )
    hits$accepted_taxon_id <- target_id
    hits$accepted_name <- NA_character_
    hits$accepted_authorship <- NA_character_
    hits$accepted_status <- NA_character_
    hits$matched_name_status <- NA_character_
    return(hits)
  }

  missing_target_ids <- setdiff(ids, as.character(accepted$taxon_id))
  if (length(missing_target_ids) > 0L) {
    warning(
      length(missing_target_ids),
      " accepted taxon ID(s) referenced by matched names were not found in ",
      "the COL SQLite database. Rebuild it from the complete NameUsage.tsv file.",
      call. = FALSE
    )
  }

  names(accepted) <- paste0("accepted_", names(accepted))
  accepted$accepted_lookup_id <- accepted$accepted_taxon_id
  hits$accepted_lookup_id <- target_id
  out <- dplyr::left_join(hits, accepted, by = "accepted_lookup_id")
  out$accepted_taxon_id <- dplyr::coalesce(
    out$accepted_taxon_id,
    out$accepted_lookup_id
  )
  out$accepted_name <- out$accepted_scientific_name
  out$accepted_authorship <- col_coalesce_columns(out, "accepted_authorship")
  out$accepted_status <- col_coalesce_columns(out, "accepted_status")
  out$matched_name_status <- out$accepted_status
  out$status <- out$accepted_status

  lineage_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  for (col in lineage_cols) {
    accepted_col <- paste0("accepted_", col)
    if (accepted_col %in% names(out)) {
      if (col == "species") {
        accepted_rank <- col_normalize_rank(col_coalesce_columns(out, "accepted_rank", "rank"))
        is_infraspecific <- accepted_rank %in% col_infraspecific_ranks()
        out[[col]] <- ifelse(
          is_infraspecific,
          dplyr::coalesce(out$accepted_scientific_name, out[[accepted_col]], out[[col]]),
          dplyr::coalesce(out[[accepted_col]], out$accepted_scientific_name, out[[col]])
        )
      } else {
        out[[col]] <- dplyr::coalesce(out[[accepted_col]], out[[col]])
      }
    }
  }

  drop_cols <- c(
    "accepted_lookup_id", "accepted_scientific_name", "accepted_rank",
    paste0("accepted_", lineage_cols)
  )
  dplyr::select(out, -dplyr::any_of(drop_cols))
}

col_drop_internal_columns <- function(x) {
  internal <- c(
    "scientific_name_norm", "scientific_name_fuzzy",
    grep("_norm$", names(x), value = TRUE)
  )
  metadata <- c(
    "alternativeid", "namealternativeid", "sourceid", "basionymid",
    "accepted_name_usage_id", "status", "accepted_status", "notho",
    "originalspelling", "uninomial"
  )
  dplyr::select(x, -dplyr::any_of(c(internal, metadata)))
}

col_coalesce_columns <- function(x, ...) {
  cols <- list(...)
  values <- lapply(cols, function(col) {
    if (col %in% names(x)) x[[col]] else rep(NA_character_, nrow(x))
  })
  do.call(dplyr::coalesce, values)
}

col_future_lapply <- function(x, fun, n_threads) {
  col_require_namespace("future")
  col_require_namespace("future.apply")

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  strategy <- if (.Platform$OS.type == "windows") future::multisession else future::multicore
  future::plan(strategy, workers = as.integer(n_threads))
  future.apply::future_lapply(x, fun, future.seed = TRUE)
}

col_group_input <- function(input, rank_cols) {
  if (length(rank_cols) == 0L) {
    return(list(input))
  }
  key <- do.call(paste, c(input[rank_cols], sep = "\r"))
  split(input, key, drop = TRUE)
}

col_group_where <- function(con, row, rank_cols) {
  clauses <- character()
  for (col in rank_cols) {
    value <- row[[col]][[1]]
    if (length(value) == 0L || is.na(value) || !nzchar(as.character(value))) {
      next
    }
    clauses <- c(
      clauses,
      sprintf(
        "%s = %s",
        DBI::dbQuoteIdentifier(con, paste0(col, "_norm")),
        DBI::dbQuoteString(con, col_normalize_name(as.character(value)))
      )
    )
  }
  clauses
}

col_sql_select <- function(con, table, selected, where) {
  sprintf(
    "SELECT %s FROM %s WHERE %s",
    paste(DBI::dbQuoteIdentifier(con, selected), collapse = ", "),
    DBI::dbQuoteIdentifier(con, table),
    paste(where[nzchar(where)], collapse = " AND ")
  )
}

col_sql_in <- function(con, column, values, quote_identifier = TRUE) {
  values <- unique(values[!is.na(values) & nzchar(as.character(values))])
  if (length(values) == 0L) {
    return("1 = 0")
  }
  lhs <- if (isTRUE(quote_identifier)) DBI::dbQuoteIdentifier(con, column) else column
  rhs <- paste(DBI::dbQuoteString(con, as.character(values)), collapse = ", ")
  sprintf("%s IN (%s)", lhs, rhs)
}

col_string_distance <- function(query, candidates, method) {
  if (requireNamespace("stringdist", quietly = TRUE)) {
    stringdist::stringdist(query, candidates, method = method)
  } else {
    as.integer(utils::adist(query, candidates))
  }
}
