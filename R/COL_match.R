#' Build a SQLite database from a Catalogue of Life TSV file
#'
#' Read a Catalogue of Life `NameUsage.tsv` file, keep only columns needed for
#' taxonomic matching, and build an indexed SQLite database with precomputed
#' lineage columns.
#'
#' @param tsv_path Path to COL `NameUsage.tsv`.
#' @param db_path Path to the SQLite database to create.
#' @param table Name of the SQLite table to write.
#' @param overwrite Logical. If `TRUE`, remove an existing database at
#'   `db_path`.
#' @param max_depth Maximum parent-chain iterations used when filling lineage
#'   columns.
#' @param chunk_size Number of rows per write chunk passed to SQLite.
#'
#' @return Invisibly returns `db_path`.
#' @export
#'
#' @examples
#' \dontrun{
#' col_build_sqlite_from_tsv("NameUsage.tsv", "col.sqlite", overwrite = TRUE)
#' }
col_build_sqlite_from_tsv <- function(tsv_path,
                                      db_path,
                                      table = "col_taxa",
                                      overwrite = FALSE,
                                      max_depth = 40,
                                      chunk_size = 100000L) {
  col_require_namespace("readr")
  col_require_namespace("dplyr")
  col_require_namespace("DBI")
  col_require_namespace("RSQLite")

  header_con <- file(tsv_path, open = "r", encoding = "UTF-8")
  header_line <- tryCatch(
    readLines(header_con, n = 1L, warn = FALSE),
    finally = close(header_con)
  )

  header <- names(readr::read_tsv(I(header_line), n_max = 0, show_col_types = FALSE))
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
      ". Use col_build_sqlite() on an in-memory table if lineage must be computed from parentID.",
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

  readr::read_tsv_chunked(
    file = tsv_path,
    callback = readr::SideEffectChunkCallback$new(callback),
    chunk_size = as.integer(chunk_size),
    col_types = keep_col_types,
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

#' Build a SQLite database from Catalogue of Life data
#'
#' Convert a Catalogue of Life taxon table to a SQLite database with normalized
#' name fields, precomputed lineage columns, and indexes for fast repeated
#' species-name lookup.
#'
#' The input table is expected to contain Catalogue of Life-style columns such as
#' `col:ID`, `col:parentID`, `col:scientificName`, `col:rank`, and `col:status`.
#' Column names may either keep the `col:` prefix or use plain names such as
#' `ID` and `scientificName`. Lineage columns are computed from `parentID` and
#' `rank`, so downstream queries can return `kingdom`, `phylum`, `class`,
#' `order`, `family`, `genus`, and `species` directly from SQLite.
#'
#' @param col_data A data frame or tibble containing COL taxon data.
#' @param db_path Path to the SQLite database to create.
#' @param table Name of the SQLite table to write.
#' @param overwrite Logical. If `TRUE`, remove an existing database at
#'   `db_path`.
#' @param max_depth Maximum parent-chain iterations used when filling lineage
#'   columns. The default is usually enough for COL ranks.
#' @param chunk_size Number of rows per write chunk passed to SQLite.
#'
#' @return Invisibly returns `db_path`.
#' @export
#'
#' @examples
#' \dontrun{
#' col_build_sqlite(col_data, "col.sqlite", overwrite = TRUE)
#' }
col_build_sqlite <- function(col_data,
                             db_path,
                             table = "col_taxa",
                             overwrite = FALSE,
                             max_depth = 40,
                             chunk_size = 100000L) {
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
#' [col_build_sqlite()] or [col_build_sqlite_from_tsv()]. Results include the
#' matched taxon record plus the full standard lineage: `kingdom`, `phylum`,
#' `class`, `order`, `family`, `genus`, and `species`.
#'
#' Exact matching is performed first using normalized scientific names. When
#' `fuzzy = TRUE`, unmatched names are compared against constrained candidate
#' sets using edit distance. The fuzzy normalization can treat common OCR or
#' typing confusions such as `i` and `l` as similar by using
#' `confusables = TRUE`.
#'
#' Constraints can be supplied in two ways. First, if `x` is a data frame and
#' contains columns named like taxonomy ranks, for example `kingdom` or
#' `family`, each row is matched only within those values. Second, fixed
#' constraints can be passed through arguments such as `kingdom = "Animalia"`.
#'
#' @param x A character vector of species names, or a data frame containing a
#'   species-name column.
#' @param db_path Path to a SQLite database created by [col_build_sqlite()].
#' @param name_col Name of the species-name column when `x` is a data frame.
#' @param table SQLite table name used by [col_build_sqlite()].
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
#' @param parallel Logical. If `TRUE`, use `doFuture` for fuzzy matching.
#' @param workers Number of parallel workers used for fuzzy matching. Defaults
#'   to 4.
#'
#' @return A tibble with original input columns prefixed by `input_` and matched
#'   COL columns. Important output columns include `match_type`,
#'   `match_distance`, `input_name_status`, `matched_name_status`, `taxon_id`,
#'   `scientific_name`, `rank`, `kingdom`, `phylum`, `class`, `order`,
#'   `family`, `genus`, `species`, `accepted_name`, and
#'   `accepted_authorship`.
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
col_match_species <- function(x,
                              db_path,
                              name_col = "species",
                              table = "col_taxa",
                              ranks = "species",
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
                              workers = 4L) {
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

  name_col <- if (is.character(x)) "species" else col_clean_names(name_col)
  input <- col_prepare_input(x, name_col)
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
  lookup_map <- input |>
    dplyr::select(dplyr::all_of(c("input_row", dedupe_cols))) |>
    dplyr::distinct(dplyr::across(dplyr::all_of(dedupe_cols)), .keep_all = TRUE) |>
    dplyr::mutate(lookup_row = dplyr::row_number())
  input <- dplyr::left_join(
    input,
    dplyr::select(lookup_map, dplyr::all_of(c(dedupe_cols, "lookup_row"))),
    by = dedupe_cols
  )
  lookup_input <- lookup_map |>
    dplyr::select(-dplyr::all_of("input_row")) |>
    dplyr::rename(input_row = "lookup_row")

  selected <- col_select_fields(fields, select)
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
      workers = workers
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

  dplyr::inner_join(input_part, out, by = "lookup_row") |>
    dplyr::arrange(.data$input_row) |>
    dplyr::select(-dplyr::all_of(c("input_row", "lookup_row"))) |>
    tibble::as_tibble()
}

col_require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required. Install it first.", call. = FALSE)
  }
}

col_required_source_columns <- function(names_clean) {
  required <- c("id", "parentid", "scientificname", "rank")
  optional <- c(
    "alternativeid", "namealternativeid", "sourceid", "basionymid", "status",
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
    specificepithet = "specific_epithet"
  )
  for (old in intersect(names(rename_map), names(taxa))) {
    names(taxa)[names(taxa) == old] <- rename_map[[old]]
  }

  taxa <- dplyr::mutate(
    taxa,
    dplyr::across(dplyr::everything(), as.character),
    rank = tolower(trimws(.data$rank)),
    scientific_name_norm = col_normalize_name(.data$scientific_name),
    scientific_name_fuzzy = col_normalize_name(.data$scientific_name, confusables = TRUE)
  )

  lineage_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  missing_lineage <- setdiff(lineage_ranks, names(taxa))
  if (length(missing_lineage) > 0L) {
    if (!isTRUE(compute_lineage)) {
      stop(
        "Chunked TSV import requires existing lineage columns: ",
        paste(lineage_ranks, collapse = ", "),
        call. = FALSE
      )
    }
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

  taxa
}

col_create_sqlite_indexes <- function(con, table) {
  lineage_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  index_cols <- intersect(
    c(
      "scientific_name_norm", "scientific_name_fuzzy", "rank",
      lineage_ranks, paste0(lineage_ranks, "_norm"), "taxon_id", "parent_id"
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
                              workers) {
  if (nrow(input) == 0L) {
    return(tibble::tibble())
  }

  rows <- lapply(seq_len(nrow(input)), function(i) input[i, , drop = FALSE])
  if (isTRUE(parallel) && length(rows) > 1L) {
    out <- col_future_lapply(rows, function(row_data) {
      worker_con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
      on.exit(DBI::dbDisconnect(worker_con), add = TRUE)
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
    }, workers = workers)
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
  where <- c(
    col_sql_in(con, "rank", ranks),
    col_sql_in(con, "substr(scientific_name_fuzzy, 1, 1)", initial, quote_identifier = FALSE),
    sprintf("length(scientific_name_fuzzy) BETWEEN %d AND %d", max(1L, min_len), max_len),
    col_group_where(con, row_data, rank_cols)
  )
  sql <- col_sql_select(con, table, selected, where)
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
  hit <- candidates[keep, , drop = FALSE]
  hit$input_row <- row_data$input_row[[1]]
  hit$match_type <- "fuzzy"
  hit$match_distance <- as.integer(raw_distances[keep])
  hit
}

col_add_accepted_names <- function(con, table, hits) {
  if (nrow(hits) == 0L || !"taxon_id" %in% names(hits)) {
    return(hits)
  }

  status <- if ("status" %in% names(hits)) hits$status else rep(NA_character_, nrow(hits))
  hits$input_name_status <- status
  target_id <- ifelse(
    !is.na(status) & status == "synonym" & !is.na(hits$parent_id) & nzchar(hits$parent_id),
    hits$parent_id,
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
  }) |>
    dplyr::bind_rows()

  if (nrow(accepted) == 0L) {
    hits$accepted_taxon_id <- target_id
    hits$accepted_name <- NA_character_
    hits$accepted_authorship <- NA_character_
    hits$accepted_status <- NA_character_
    return(hits)
  }

  names(accepted) <- paste0("accepted_", names(accepted))
  accepted$accepted_lookup_id <- accepted$accepted_taxon_id
  hits$accepted_lookup_id <- target_id
  out <- dplyr::left_join(hits, accepted, by = "accepted_lookup_id")
  out$accepted_taxon_id <- dplyr::coalesce(out$accepted_taxon_id, out$taxon_id)
  out$accepted_name <- dplyr::coalesce(out$accepted_scientific_name, out$scientific_name)
  out$accepted_authorship <- col_coalesce_columns(out, "accepted_authorship", "authorship")
  out$accepted_status <- col_coalesce_columns(out, "accepted_status", "status")
  out$matched_name_status <- out$accepted_status
  out$status <- out$accepted_status

  lineage_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  for (col in lineage_cols) {
    accepted_col <- paste0("accepted_", col)
    if (accepted_col %in% names(out)) {
      if (col == "species") {
        out[[col]] <- dplyr::coalesce(out[[accepted_col]], out$accepted_scientific_name, out[[col]])
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
    "status", "accepted_status", "notho",
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

col_future_lapply <- function(x, fun, workers) {
  col_require_namespace("doFuture")
  col_require_namespace("future")
  col_require_namespace("foreach")

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  strategy <- if (.Platform$OS.type == "windows") future::multisession else future::multicore
  future::plan(strategy, workers = as.integer(workers))
  doFuture::registerDoFuture()

  foreach::`%dopar%`(
    foreach::foreach(
      i = x,
      .export = c(
        "col_fuzzy_match_one", "col_string_distance", "col_sql_in",
        "col_sql_select", "col_group_where", "col_normalize_name"
      ),
      .packages = character()
    ),
    fun(i)
  )
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
