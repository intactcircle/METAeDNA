# NCBI taxonomy database construction and querying
.validate_accession_index_columns <- function(accession_index_columns) {
  valid_accession_indexes <- c("accession", "accession.version", "taxid")
  if (!is.character(accession_index_columns) ||
      anyNA(accession_index_columns) ||
      any(!accession_index_columns %in% valid_accession_indexes)) {
    stop(
      "accession_index_columns must contain only: ",
      paste(valid_accession_indexes, collapse = ", "),
      call. = FALSE
    )
  }
  unique(accession_index_columns)
}

.create_taxonomy_indexes <- function(con) {
  existing_tables <- DBI::dbListTables(con)
  missing_tables <- setdiff(c("nodes", "names", "merged"), existing_tables)
  if (length(missing_tables) > 0L) {
    stop(
      "Cannot create taxonomy indexes; missing table(s): ",
      paste(missing_tables, collapse = ", "),
      call. = FALSE
    )
  }

  show_time("Indexing [taxonomy] table ...")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_nodes_parent ON nodes(parent_tax_id);")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_names_tax ON names(tax_id);")
  DBI::dbExecute(con, "CREATE UNIQUE INDEX IF NOT EXISTS idx_merged_old_tax ON merged(old_tax_id);")
  if ("name_aliases" %in% existing_tables) {
    DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_name_aliases_name ON name_aliases(name);")
    DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_name_aliases_tax ON name_aliases(tax_id);")
  }
  if ("species_equivalents" %in% existing_tables) {
    DBI::dbExecute(
      con,
      paste(
        "CREATE UNIQUE INDEX IF NOT EXISTS idx_species_equivalents_source",
        "ON species_equivalents(source_tax_id);"
      )
    )
  }
}

.ncbi_canonical_species_name <- function(x) {
  x <- trimws(as.character(x))
  qualifier <- paste0(
    "\\s+(?:s\\s*\\.?\\s*l\\s*\\.?|sensu\\s+lato)",
    "(?:\\s+.*)?$"
  )
  trimws(sub(qualifier, "", x, ignore.case = TRUE, perl = TRUE))
}

.ncbi_unbracket_name <- function(x) {
  gsub("\\[([^][]+)\\]", "\\1", as.character(x), perl = TRUE)
}

.resolve_ncbi_taxids <- function(con, taxids) {
  resolved <- suppressWarnings(as.integer(taxids))
  if (length(resolved) == 0L) {
    return(resolved)
  }

  redirect_with_table <- function(values, table, old_col, new_col) {
    query_ids <- unique(values[!is.na(values)])
    if (length(query_ids) == 0L || !table %in% DBI::dbListTables(con)) {
      return(values)
    }
    redirects <- DBI::dbGetQuery(
      con,
      paste0(
        "SELECT ", old_col, ", ", new_col, " FROM ", table,
        " WHERE ", old_col, " IN (",
        paste(DBI::dbQuoteLiteral(con, query_ids), collapse = ","),
        ")"
      )
    )
    redirect_match <- match(values, redirects[[old_col]])
    replace <- !is.na(redirect_match)
    values[replace] <- redirects[[new_col]][redirect_match[replace]]
    values
  }

  resolved <- redirect_with_table(
    resolved, "merged", "old_tax_id", "new_tax_id"
  )
  redirect_with_table(
    resolved,
    "species_equivalents",
    "source_tax_id",
    "canonical_tax_id"
  )
}

.create_accession_indexes <- function(con, accession_index_columns = c("accession", "accession.version")) {
  accession_index_columns <- .validate_accession_index_columns(accession_index_columns)
  if (length(accession_index_columns) == 0L) {
    return(invisible(NULL))
  }

  if (!"accession2taxid" %in% DBI::dbListTables(con)) {
    stop(
      "Cannot create accession indexes; missing table: accession2taxid.",
      call. = FALSE
    )
  }

  show_time("Indexing [accession2taxid] table ...")
  accession_index_sql <- c(
    accession = "CREATE INDEX IF NOT EXISTS idx_accession2taxid_accession ON accession2taxid(accession);",
    `accession.version` = "CREATE INDEX IF NOT EXISTS idx_accession2taxid_accession_version ON accession2taxid(`accession.version`);",
    taxid = "CREATE INDEX IF NOT EXISTS idx_accession2taxid_taxid ON accession2taxid(taxid);"
  )
  for (column in accession_index_columns) {
    DBI::dbExecute(con, accession_index_sql[[column]])
  }
  invisible(NULL)
}

#' Build indexes for an existing NCBI taxonomy SQLite database
#'
#' Add indexes after the database tables have already been imported. This is
#' useful when \code{build_ncbi_taxonomy_database(create_index = FALSE)} was used to finish the
#' large accession/WGS import first and postpone the slower indexing step.
#'
#' @param sql_path (character) path to an existing SQLite database created by
#'   \code{build_ncbi_taxonomy_database()}.
#' @param taxonomy (logical) if \code{TRUE}, create indexes for \code{nodes},
#'   \code{names}, \code{name_aliases}, \code{species_equivalents}, and
#'   \code{merged}.
#' @param accession (logical) if \code{TRUE}, create indexes for
#'   \code{accession2taxid}.
#' @param accession_index_columns (character vector) accession-table columns to
#'   index. The default indexes \code{accession} and
#'   \code{accession.version}; add \code{"taxid"} only if reverse taxid lookup
#'   is required.
#'
#' @return Invisibly returns \code{NULL}.
#' @export
build_database_index <- function(sql_path,
                                 taxonomy = FALSE,
                                 accession = TRUE,
                                 accession_index_columns = c("accession", "accession.version")) {
  if (!is.character(sql_path) || length(sql_path) != 1L || is.na(sql_path)) {
    stop("sql_path must be a single non-NA character value.", call. = FALSE)
  }
  if (!file.exists(sql_path)) {
    stop("SQLite database does not exist: ", sql_path, call. = FALSE)
  }
  accession_index_columns <- .validate_accession_index_columns(accession_index_columns)

  start_time <- Sys.time()
  on.exit(show_time(secs = as.numeric(Sys.time() - start_time, units = "secs"), mode = 1))

  con <- DBI::dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  if (isTRUE(taxonomy)) {
    .create_taxonomy_indexes(con)
  }
  if (isTRUE(accession)) {
    .create_accession_indexes(con, accession_index_columns)
  }
  show_time("Done.", mode = 0)
  invisible(NULL)
}

#' Build an NCBI taxonomy and accession database
#'
#' This function builds an SQLite database from NCBI taxonomy resources. It can
#' optionally download and parse the NCBI taxdump archive (nodes, names, and
#' merged tax IDs) and
#' the large \code{accession2taxid} mapping files, then write them into an
#' SQLite database for fast querying by downstream functions such as
#' \code{get_taxid()} and \code{get_taxonomy()}. The function supports
#' incremental construction of taxonomy-only, accession-only, or combined
#' databases, and can overwrite an existing database if requested.
#'
#' @param taxonomy (logical) if \code{TRUE}, create taxonomy tables from
#'   \code{nodes.dmp}, \code{names.dmp}, and \code{merged.dmp} in the NCBI
#'   taxdump archive. These tables are written as \code{nodes}, \code{names},
#'   \code{merged}, \code{name_aliases}, and \code{species_equivalents} in the
#'   SQLite database.
#'   Alternatively, the files can be downloaded directly from NCBI at: \code{https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/}, and put this file to
#'   your working space (specified in \code{taxdump_dir}).
#'
#' @param accession (logical) if \code{TRUE}, create an \code{accession2taxid}
#'   table using the \code{nucl_gb.accession2taxid.gz} file provided by NCBI.
#'
#' @param include_WGS (logical) if \code{TRUE}, also include WGS accessions by
#'   importing \code{nucl_wgs.accession2taxid.gz} and appending them into the
#'   \code{accession2taxid} table. This substantially increases database size and processing time.
#'
#' @param taxdump_dir (character) cache directory for all downloaded source
#'   files. Existing, valid files in this directory are reused; missing or
#'   invalid required files are downloaded when \code{download = TRUE}. The
#'   directory must have sufficient free disk space (typically at least
#'   10-120 GB, depending on whether accession tables and WGS data are included).
#'
#' @param merged_path (character or \code{NULL}) optional path to an existing
#'   \code{merged.dmp} file. The default \code{NULL} automatically discovers
#'   \code{merged.dmp} in \code{taxdump_dir}; users normally do not need to set
#'   this argument. An explicit path may be supplied to use a shared file, such
#'   as \code{"reference_databases/ncbi_taxonomy/merged.dmp"},
#'   without copying it into \code{taxdump_dir}.
#'
#' @param sql_path (character) file path where the SQLite database will be
#'   created. Using an absolute path is recommended. The target file system
#'   should provide at least 120 GB of free space when building a full
#'   accession2taxid database with WGS data and indexed.
#'
#' @param overwrite (logical) whether existing data may be replaced. If
#'   \code{FALSE} and \code{sql_path} exists, no changes are made.
#'
#' @param overwrite_scope (character) controls what \code{overwrite = TRUE}
#'   replaces when the database already exists. \code{"database"} removes the
#'   entire SQLite file before rebuilding, preserving the historical behavior.
#'   \code{"selected"} replaces only the components selected by
#'   \code{taxonomy} and \code{accession}; unselected tables are preserved.
#'   Thus, setting \code{taxonomy = TRUE}, \code{accession = FALSE}, and
#'   \code{overwrite_scope = "selected"} updates taxonomy without rebuilding
#'   the large \code{accession2taxid} table.
#'
#' @param download (logical) if \code{TRUE}, download required files that are
#'   missing or whose local size differs from the current NCBI
#'   \code{Content-Length} by more than 5 percent. Downloads are written to
#'   temporary \code{.part} files and moved to their NCBI file names only after
#'   the download succeeds and passes the same size comparison. Large gzip
#'   files are not decompressed solely for preliminary validation. Their
#'   complete compressed streams and expected contents are checked during the
#'   subsequent extraction or C++ import, so each large input is decompressed
#'   only once. If the remote size cannot be obtained, an existing nonempty
#'   local file is retained with a warning.
#'
#' @param resume (logical) if \code{TRUE}, resume an interrupted download from
#'   its \code{.part} file using \code{curl::multi_download()}.
#'
#' @param create_index (logical) if \code{TRUE}, create indexes on key columns
#'   (e.g., \code{parent_tax_id} in \code{nodes}, \code{tax_id} in \code{names},
#'   \code{old_tax_id} in \code{merged}, names in \code{name_aliases}, and
#'   \code{source_tax_id} in \code{species_equivalents}, as well as
#'   \code{accession} or \code{taxid} in \code{accession2taxid}) to improve
#'   query performance at the cost of additional disk space and initial build
#'   time.
#'
#' @param chunk_size (positive integer) maximum number of accession records
#'   retained by the C++ gzip parser before one batch is passed to RSQLite. The
#'   default, \code{500000}, balances speed and memory use on most systems.
#'   Reduce it to lower peak memory use or increase it to reduce append
#'   overhead. The complete accession file is never loaded into R memory.
#'
#' @param accession_index_columns (character vector) accession-table columns
#'   to index when \code{create_index = TRUE}. The default indexes
#'   \code{accession} and \code{accession.version}, which are used by
#'   \code{get_taxid()}. Add \code{"taxid"} only when reverse taxid lookups are
#'   required; omitting that very large third index substantially reduces WGS
#'   build time and disk use. Use \code{character(0)} to skip accession indexes
#'   while still creating taxonomy indexes.
#'
#' @details
#'
#' Building the accession-to-taxid table, especially with WGS data, requires
#' substantial disk space and I/O time. Progress messages and approximate
#' record counts are printed during import. Setting \code{include_WGS = TRUE}
#' implies \code{accession = TRUE}; both \code{nucl_gb.accession2taxid.gz} and
#' \code{nucl_wgs.accession2taxid.gz} are then required and imported. Accession
#' files are read directly from gzip in chunks using the fixed NCBI four-column
#' schema, including \code{gi}. Records are first written to a staging table
#' inside a transaction, and the existing table is replaced only after every
#' requested file has been imported successfully.
#' Because the staging table is disposable until activation, its bulk import
#' uses SQLite's no-journal mode for faster writes. The normal journal mode is
#' restored before the staging table is atomically activated.
#'
#' With \code{overwrite_scope = "selected"}, taxonomy tables are also built as
#' staging tables and swapped into place together. A failed taxonomy update
#' therefore does not remove the existing taxonomy or accession tables.
#'
#' Species-level NCBI records qualified with \code{s.l.} or
#' \code{sensu lato} may represent separately sequenced samples of an existing
#' binomial, for example \code{Gnetum hainanense s.l. Won 511}. When the
#' unqualified binomial identifies exactly one species node, the taxonomy build
#' records the qualified tax ID in \code{species_equivalents}. Subsequent
#' accession, tax-ID, and taxon-name queries resolve it to that canonical
#' species. Named infraspecific taxa such as varieties and subspecies are not
#' changed by this rule.
#'
#' NCBI may enclose a provisional genus in square brackets, for example
#' \code{[Polygonum] chinense var. procumbens}. If the same tax ID also contains
#' the exactly corresponding unbracketed synonym, that verified form is used as
#' the displayed scientific name. Its original tax ID and infraspecific rank
#' are preserved, and both forms remain searchable in \code{name_aliases}.
#'
#' @return
#' Invisibly returns \code{NULL}. The main result is the SQLite database file
#' created at \code{sql_path}, containing some or all of the following tables:
#' \code{nodes}, \code{names}, \code{name_aliases},
#' \code{species_equivalents}, \code{merged}, and \code{accession2taxid},
#' depending on the arguments supplied.
#' @importFrom DBI dbAppendTable dbBegin dbCommit dbConnect dbDisconnect
#'   dbExecute dbRollback dbWriteTable
#' @importFrom RSQLite SQLite
#' @importFrom readr col_character col_integer cols read_delim read_tsv
#' @importFrom curl curl_download
#' @importFrom utils tail untar
#' @importFrom stringr str_remove
#'
#' @examples
#' \dontrun{
#' # Build taxonomy-only database in the current directory
#' build_ncbi_taxonomy_database(
#'   taxonomy   = TRUE,
#'   accession  = FALSE,
#'   include_WGS = FALSE,
#'   taxdump_dir = ".",
#'   sql_path    = "./taxdump.sqlite",
#'   overwrite   = FALSE
#' )
#'
#' # Build taxonomy + accession2taxid (GenBank) database with overwrite
#' build_ncbi_taxonomy_database(
#'   taxonomy   = TRUE,
#'   accession  = TRUE,
#'   include_WGS = FALSE,
#'   taxdump_dir = "reference_databases/ncbi_taxonomy",
#'   sql_path    = "reference_databases/ncbi_taxonomy/taxonomy.sqlite",
#'   overwrite   = TRUE
#' )
#'
#' # Continue the workflow: accession -> taxid -> taxonomy
#' db <- "reference_databases/ncbi_taxonomy/taxonomy.sqlite"
#' accession_taxids <- get_taxid(
#'   sql_path = db,
#'   accessions = c("NC_045512.2", "NC_000001.11")
#' )
#' accession_lineages <- get_taxonomy(
#'   sql_path = db,
#'   tax_ids = accession_taxids$taxid
#' )
#'
#' # Refresh taxonomy and redirects without rebuilding accession2taxid
#' build_ncbi_taxonomy_database(
#'   taxonomy = TRUE,
#'   accession = FALSE,
#'   taxdump_dir = "reference_databases/ncbi_taxonomy",
#'   sql_path = db,
#'   overwrite = TRUE,
#'   overwrite_scope = "selected"
#' )
#' }
#'
#' @seealso
#' \code{\link{get_taxid}} for mapping accessions to NCBI tax IDs, and
#' \code{\link{get_taxonomy}} for retrieving taxonomic lineages from the
#' database created by this function.
#'
#' @export
build_ncbi_taxonomy_database <- function(taxonomy = TRUE, accession = FALSE, include_WGS = FALSE, taxdump_dir = ".", sql_path = "./taxdump.sqlite", overwrite = F, download = T,resume = T,create_index = T, chunk_size = 1000000L, merged_path = NULL, overwrite_scope = c("database", "selected"), accession_index_columns = c("accession", "accession.version")){
  start_time <- Sys.time()
  on.exit(show_time(secs = as.numeric(Sys.time() - start_time, units = "secs"), mode = 1))
  overwrite_scope <- match.arg(overwrite_scope)
  accession_index_columns <- .validate_accession_index_columns(accession_index_columns)
  if (!is.null(merged_path) &&
      (!is.character(merged_path) || length(merged_path) != 1L ||
       is.na(merged_path) || !nzchar(merged_path))) {
    stop("merged_path must be NULL or a single non-empty file path.", call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      is.na(chunk_size) || !is.finite(chunk_size) ||
      chunk_size < 1 || chunk_size > .Machine$integer.max ||
      chunk_size != floor(chunk_size)) {
    stop("chunk_size must be a single positive integer.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)
  if (!dir.exists(taxdump_dir) && !dir.create(taxdump_dir, recursive = TRUE)) {
    stop("Could not create taxdump_dir: ", taxdump_dir)
  }

  local_file_size <- function(path) {
    info <- file.info(path)
    if (!file.exists(path) || nrow(info) != 1L || is.na(info$size) ||
        info$size < 1) {
      return(NA_real_)
    }
    as.numeric(info$size)
  }

  remote_file_size <- function(url) {
    tryCatch({
      handle <- curl::new_handle(
        nobody = TRUE,
        followlocation = TRUE,
        accept_encoding = "identity",
        connecttimeout = 15,
        timeout = 30
      )
      response <- curl::curl_fetch_memory(url, handle = handle)
      if (response$status_code < 200L || response$status_code >= 400L) {
        return(NA_real_)
      }
      headers <- curl::parse_headers_list(response$headers)
      size <- suppressWarnings(as.numeric(headers[["content-length"]]))
      if (length(size) != 1L || is.na(size) || !is.finite(size) || size < 1) {
        return(NA_real_)
      }
      size
    }, error = function(e) NA_real_)
  }

  validate_remote_size <- function(path, url, tolerance = 0.05) {
    local_size <- local_file_size(path)
    if (is.na(local_size)) return(FALSE)

    remote_size <- remote_file_size(url)
    if (is.na(remote_size)) {
      warning(
        "Could not obtain the remote file size; using the existing nonempty file: ",
        path,
        call. = FALSE
      )
      return(TRUE)
    }

    relative_difference <- abs(local_size - remote_size) / remote_size
    is.finite(relative_difference) && relative_difference <= tolerance
  }

  validate_taxdump <- function(path, url) {
    validate_remote_size(path, url)
  }

  validate_dmp <- function(path) {
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
      return(FALSE)
    }
    first <- tryCatch(readLines(path, n = 1L, warn = FALSE),
                      error = function(e) character())
    tail_lines <- tryCatch({
      con <- file(path, open = "rb")
      on.exit(close(con))
      seek(con, where = max(0, file.info(path)$size - 8192), origin = "start")
      readLines(con, warn = FALSE)
    }, error = function(e) character())
    last <- tail(tail_lines[nzchar(tail_lines)], 1L)
    length(first) == 1L && length(last) == 1L &&
      grepl("\t\\|", first) && grepl("\t\\|$", last)
  }

  validate_accession_gzip <- function(path, url) {
    validate_remote_size(path, url)
  }

  get_required_file <- function(filename, url, validator) {
    destination <- file.path(taxdump_dir, filename)
    if (validator(destination, url)) {
      show_time(sprintf("Use existing [%s].", filename))
      return(destination)
    }
    if (!download) {
      stop(
        "Required file is missing, empty, or differs from the remote size by ",
        "more than 5%: ",
        destination
      )
    }

    part <- paste0(destination, ".part")
    if (file.exists(destination)) {
      warning("Ignoring invalid file: ", destination, call. = FALSE)
    }
    show_time(sprintf("Download [%s]...", filename))
    curl::multi_download(url, destfiles = part, resume = resume)
    if (!validator(part, url)) {
      stop(
        "Downloaded file size differs from the remote file by more than 5%; ",
        "partial file retained for resuming: ",
        part
      )
    }
    if (file.exists(destination) && unlink(destination) != 0L) {
      stop("Could not replace invalid file: ", destination)
    }
    if (!file.rename(part, destination)) {
      stop("Downloaded file could not be moved into place: ", destination)
    }
    destination
  }

  accession <- isTRUE(accession) || isTRUE(include_WGS)
  show_time("Building the SQLite database...")
  if(file.exists(sql_path) && !overwrite){
    warning(
      "The database already exists. Set overwrite = TRUE and choose ",
      "overwrite_scope to update it.",
      call. = FALSE
    )
    return(invisible(NULL))
  }
  if(file.exists(sql_path) && overwrite && overwrite_scope == "database"){
    if (unlink(sql_path) != 0L) {
      stop("Could not remove the existing database: ", sql_path)
    }
  } else if (file.exists(sql_path) && overwrite) {
    show_time("Updating only the selected database components; unselected tables will be preserved.")
  }
  if(taxonomy){
    nodes_path <- file.path(taxdump_dir, "nodes.dmp")
    names_path <- file.path(taxdump_dir, "names.dmp")
    merged_path_resolved <- if (is.null(merged_path)) {
      file.path(taxdump_dir, "merged.dmp")
    } else {
      path.expand(merged_path)
    }
    cached_merged_required <- is.null(merged_path)
    if (!validate_dmp(nodes_path) || !validate_dmp(names_path) ||
        (cached_merged_required && !validate_dmp(merged_path_resolved))) {
      taxdump_path <- get_required_file(
        "taxdump.tar.gz",
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        validate_taxdump
      )
      utils::untar(taxdump_path, exdir = taxdump_dir)
    }
    if (!validate_dmp(nodes_path) || !validate_dmp(names_path)) {
      stop("taxdump.tar.gz did not produce valid nodes.dmp and names.dmp files.")
    }
    if (!validate_dmp(merged_path_resolved)) {
      stop("merged.dmp is missing or invalid: ", merged_path_resolved)
    }
    show_time(sprintf(
      "Use [nodes.dmp], [names.dmp], and merged redirects from [%s].",
      merged_path_resolved
    ))

    # nodes dmp structure ------------------------------
    # tax_id					-- node id in GenBank taxonomy database
    # parent tax_id				-- parent node id in GenBank taxonomy database
    # rank					-- rank of this node (domain, kingdom, ...)
    # embl code				-- locus-name prefix; not unique
    # division id				-- see division.dmp file
    # inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
    # genetic code id				-- see gencode.dmp file
    # inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
    # mitochondrial genetic code id		-- see gencode.dmp file
    # inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
    # GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
    # hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
    # comments				-- free-text comments and citations
    # ----------------------------------------------------
    # 1) read nodes ---------------------------------------
    nodes_cols <- c(
      "tax_id", "parent_tax_id", "rank", "embl_code", "division_id",
      "inherited_div_flag", "genetic_code_id", "inherited_gc_flag",
      "mito_genetic_code_id", "inherited_mgc_flag", "genbank_hidden_flag",
      "subtree_root_flag", "comments"
    )
    nodes <- read_delim(file.path(taxdump_dir, "nodes.dmp"), delim = "\t|\t", col_names = nodes_cols, show_col_types = F, progress = F) %>%
      mutate(
        comments = comments %>% str_remove("\t[|]$")
      )
    # names dmp structure ------------------------
    # tax_id					-- the id of node associated with this name
    # name_txt				-- name itself
    # unique name				-- the unique variant of this name if name not unique
    # name class				-- (synonym, common name, ...)
    # ----------------------------------
    # 2) read names.dmp -------------------------------------------
    names_cols <- c("tax_id", "name_txt", "unique_name", "name_class")
    names_raw <- read_delim(
      file.path(taxdump_dir, "names.dmp"),
      delim = "\t|\t",
      col_names = names_cols,
      show_col_types = FALSE,
      progress = FALSE
    ) %>%
      mutate(name_class = name_class %>% str_remove("\t[|]$"))
    names <- names_raw %>%
      filter(name_class == "scientific name") %>%
      select(tax_id = tax_id, name = name_txt)
    name_aliases <- names_raw %>%
      filter(
        name_class %in% c(
          "scientific name", "synonym", "equivalent name",
          "genbank synonym", "anamorph", "teleomorph"
        )
      ) %>%
      select(tax_id = tax_id, name = name_txt, name_class = name_class) %>%
      distinct(tax_id, name, .keep_all = TRUE)

    # NCBI sometimes brackets a provisional genus in the scientific name while
    # providing the same unbracketed name as a synonym for the same tax ID.
    # Prefer that verified display form without changing the tax ID or rank.
    unbracketed_aliases <- name_aliases %>%
      filter(
        name_class %in% c(
          "synonym", "equivalent name", "genbank synonym"
        )
      )
    unbracketed_alias_keys <- paste(
      unbracketed_aliases$tax_id,
      unbracketed_aliases$name,
      sep = "\r"
    )
    bracket_replacements <- names %>%
      filter(grepl("\\[[^]]+\\]", name)) %>%
      mutate(display_name = .ncbi_unbracket_name(name)) %>%
      filter(
        paste(tax_id, display_name, sep = "\r") %in%
          unbracketed_alias_keys
      ) %>%
      distinct(tax_id, .keep_all = TRUE) %>%
      select(tax_id, display_name)
    if (nrow(bracket_replacements) > 0L) {
      replacement_match <- match(names$tax_id, bracket_replacements$tax_id)
      replace <- !is.na(replacement_match)
      names$name[replace] <-
        bracket_replacements$display_name[replacement_match[replace]]
    }
    rm(names_raw)

    # Map NCBI species nodes qualified as sensu lato to an unambiguous
    # canonical species node with the same binomial.
    species_names <- names %>%
      inner_join(
        nodes %>% select(tax_id, rank),
        by = "tax_id"
      ) %>%
      filter(rank == "species") %>%
      mutate(canonical_name = .ncbi_canonical_species_name(name))
    canonical_targets <- species_names %>%
      filter(name == canonical_name) %>%
      group_by(canonical_name) %>%
      filter(dplyr::n_distinct(tax_id) == 1L) %>%
      ungroup() %>%
      select(canonical_name, canonical_tax_id = tax_id)
    species_equivalents <- species_names %>%
      filter(name != canonical_name) %>%
      inner_join(canonical_targets, by = "canonical_name") %>%
      filter(tax_id != canonical_tax_id) %>%
      dplyr::transmute(
        source_tax_id = tax_id,
        canonical_tax_id = canonical_tax_id,
        source_name = name,
        canonical_name = canonical_name
      ) %>%
      distinct(source_tax_id, .keep_all = TRUE)
    rm(species_names, canonical_targets)

    # 3) read redirects from obsolete tax IDs to their current tax IDs -------
    merged <- read_delim(
      merged_path_resolved,
      delim = "\t|\t",
      col_names = c("old_tax_id", "new_tax_id"),
      col_types = readr::cols(
        old_tax_id = readr::col_integer(),
        new_tax_id = readr::col_character()
      ),
      show_col_types = FALSE,
      progress = FALSE
    ) %>%
      mutate(new_tax_id = as.integer(str_remove(new_tax_id, "\t[|]$")))

    # 4) stage and atomically replace taxonomy tables ---------------
    con <- dbConnect(SQLite(), sql_path)
    on.exit(if (DBI::dbIsValid(con)) dbDisconnect(con), add = TRUE)
    dbExecute(con, "PRAGMA journal_mode=DELETE;")

    taxonomy_staging <- c(
      nodes = "nodes_staging",
      names = "names_staging",
      name_aliases = "name_aliases_staging",
      species_equivalents = "species_equivalents_staging",
      merged = "merged_staging"
    )
    for (table in unname(taxonomy_staging)) {
      dbExecute(con, sprintf("DROP TABLE IF EXISTS %s;", table))
    }
    staging_error <- tryCatch({
      dbWriteTable(con, taxonomy_staging[["nodes"]], nodes, overwrite = TRUE)
      dbWriteTable(con, taxonomy_staging[["names"]], names, overwrite = TRUE)
      dbWriteTable(con, taxonomy_staging[["name_aliases"]], name_aliases, overwrite = TRUE)
      dbWriteTable(
        con,
        taxonomy_staging[["species_equivalents"]],
        species_equivalents,
        overwrite = TRUE
      )
      dbWriteTable(con, taxonomy_staging[["merged"]], merged, overwrite = TRUE)
      NULL
    }, error = function(e) e)
    if (!is.null(staging_error)) {
      stop(
        "Failed to stage taxonomy data; existing database tables were preserved: ",
        conditionMessage(staging_error),
        call. = FALSE
      )
    }

    DBI::dbBegin(con)
    taxonomy_swap_error <- tryCatch({
      for (table in names(taxonomy_staging)) {
        dbExecute(con, sprintf("DROP TABLE IF EXISTS %s;", table))
        dbExecute(
          con,
          sprintf(
            "ALTER TABLE %s RENAME TO %s;",
            taxonomy_staging[[table]],
            table
          )
        )
      }
      DBI::dbCommit(con)
      NULL
    }, error = function(e) {
      try(DBI::dbRollback(con), silent = TRUE)
      e
    })
    if (!is.null(taxonomy_swap_error)) {
      stop(
        "Failed to activate staged taxonomy data; existing tables were preserved: ",
        conditionMessage(taxonomy_swap_error),
        call. = FALSE
      )
    }

    if(create_index){
      .create_taxonomy_indexes(con)
    }
    dbDisconnect(con)
  }
  if(accession){
    warning(
      "Building accession2taxid may require more than 120 GB of free disk space.",
      call. = FALSE,
      immediate. = TRUE
    )
    message(
      "Existing nucl_gb.accession2taxid.gz files in taxdump_dir are reused."
    )
    if(include_WGS){
      message(
        "include_WGS = TRUE also requires nucl_wgs.accession2taxid.gz."
      )
    }
    get_required_file(
      "nucl_gb.accession2taxid.gz",
      "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
      validate_accession_gzip
    )
  }
  if(include_WGS){
    get_required_file(
      "nucl_wgs.accession2taxid.gz",
      "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
      validate_accession_gzip
    )
  }
  if(accession){
    show_time("Constructing the accession2id table...")
    con <- dbConnect(SQLite(), sql_path)
    on.exit(if (DBI::dbIsValid(con)) dbDisconnect(con), add = TRUE)

    dbExecute(con, "PRAGMA journal_mode=OFF;")
    dbExecute(con, "PRAGMA synchronous=OFF;")
    dbExecute(con, "PRAGMA locking_mode=EXCLUSIVE;")
    dbExecute(con, "PRAGMA temp_store=MEMORY;")
    dbExecute(con, "PRAGMA cache_size=-65536;")

    staging_table <- "accession2taxid_staging"
    dbExecute(con, sprintf("DROP TABLE IF EXISTS %s;", staging_table))
    dbExecute(
      con,
      sprintf(
        paste(
          "CREATE TABLE %s (",
          "accession TEXT NOT NULL,",
          "`accession.version` TEXT NOT NULL,",
          "taxid INTEGER NOT NULL,",
          "gi INTEGER NOT NULL",
          ");"
        ),
        staging_table
      )
    )

    accession_files <- file.path(taxdump_dir, "nucl_gb.accession2taxid.gz")
    if(include_WGS){
      accession_files <- c(
        accession_files,
        file.path(taxdump_dir, "nucl_wgs.accession2taxid.gz")
      )
    }
    n_processed <- 0

    import_error <- tryCatch({
      DBI::dbBegin(con)
      for(accession_file in accession_files){
        file_records <- 0
        callback <- function(chunk){
          DBI::dbAppendTable(con, staging_table, chunk)
          file_records <<- file_records + nrow(chunk)
          n_processed <<- n_processed + nrow(chunk)
          show_time(
            sprintf(
              "Processing records in [%s]: %d (total: %d)...\r",
              basename(accession_file), file_records, n_processed
            ),
            mode = 2
          )
          invisible(NULL)
        }
        stream_accession2taxid_gz(
          input_file = accession_file,
          chunk_size = chunk_size,
          callback = callback
        )
      }
      DBI::dbCommit(con)
      NULL
    }, error = function(e){
      try(DBI::dbRollback(con), silent = TRUE)
      e
    })

    # Previous implementation retained for reference. It parsed each gzip file
    # into R tibbles before appending them to SQLite, which added substantial
    # conversion overhead for WGS data:
    #
    # callback <- readr::DataFrameCallback$new(function(chunk, pos) {
    #   DBI::dbAppendTable(con, staging_table, chunk)
    # })
    # readr::read_tsv_chunked(
    #   file = accession_file,
    #   callback = callback,
    #   chunk_size = chunk_size,
    #   col_types = accession_col_types,
    #   progress = FALSE,
    #   show_col_types = FALSE
    # )
    cat("\n")

    if(!is.null(import_error)){
      dbExecute(con, sprintf("DROP TABLE IF EXISTS %s;", staging_table))
      dbExecute(con, "PRAGMA journal_mode=DELETE;")
      stop(
        "Failed to import accession2taxid data; the existing table was preserved: ",
        conditionMessage(import_error),
        call. = FALSE
      )
    }
    if(n_processed == 0){
      dbExecute(con, sprintf("DROP TABLE IF EXISTS %s;", staging_table))
      dbExecute(con, "PRAGMA journal_mode=DELETE;")
      stop("No accession2taxid records were imported.", call. = FALSE)
    }

    dbExecute(con, "PRAGMA journal_mode=DELETE;")
    DBI::dbBegin(con)
    swap_error <- tryCatch({
      dbExecute(con, "DROP TABLE IF EXISTS accession2taxid;")
      dbExecute(
        con,
        sprintf("ALTER TABLE %s RENAME TO accession2taxid;", staging_table)
      )
      DBI::dbCommit(con)
      NULL
    }, error = function(e){
      try(DBI::dbRollback(con), silent = TRUE)
      e
    })
    if(!is.null(swap_error)){
      stop(
        "Accession data were imported but the staging table could not be activated: ",
        conditionMessage(swap_error),
        call. = FALSE
      )
    }
    dbDisconnect(con)
  }

  if(accession && create_index && length(accession_index_columns) > 0L){
    con <- dbConnect(SQLite(), sql_path)
    .create_accession_indexes(con, accession_index_columns)
    dbDisconnect(con)
  }

  #
  show_time("Done.", mode = 0)
}

# Query and retrieve taxonomic information-------------------------------------
#' Retrieve Taxonomic Lineages from a Local NCBI Taxonomy SQLite Database
#'
#' This function queries a local SQLite database containing NCBI taxonomy
#' tables (`nodes` and `names`) and returns taxonomic lineages for a set of
#' tax IDs or taxon names. The function performs a recursive SQL query to
#' traverse parent nodes until reaching the root, then extracts and reshapes
#' lineage information for user-defined taxonomic ranks.
#'
#' @param sql_path (character) File path to the SQLite database containing
#'   NCBI taxonomy tables (`nodes`, `names`).
#' @param tax_ids (integer vector) NCBI taxonomy IDs to query.
#'   Cannot be used together with `tax_names`.
#' @param tax_names (character vector). Scientific names to query.
#'   The function will internally retrieve the corresponding tax IDs.
#' @param dereplicate (logical) If `TRUE`, duplicated query entries are removed
#'   after lineage construction. Default is `TRUE`.
#' @param rank_name (character vector) Ordered taxonomic ranks to extract.
#'   Must be given from highest to lowest rank. The default is
#'   `c("kingdom","phylum","class","order","family","genus","species")`.
#' @param show_elapsed (logical) If `TRUE`, elapsed time is printed using
#'   the internal `show_time()` function. Default is `TRUE`.
#' @param keep_input_taxid (logical) if \code{TRUE}, include an
#'   \code{input_taxid} column containing the originally requested IDs.
#'   \code{query_taxid} always contains the current ID after applying
#'   \code{merged.dmp} redirects. Default is \code{FALSE}.
#'
#' @details
#' The function constructs a recursive SQL query using a Common Table
#' Expression (CTE) to follow taxonomic lineages upward from each query tax ID.
#' Returned results include only user-defined ranks. When `tax_names` is used,
#' the function also returns the matched tax IDs for each supplied name.
#' Common infraspecific NCBI ranks, including \code{subspecies},
#' \code{varietas}, cultivar groups, strains, genotypes, and unranked
#' subgroups, are returned in the \code{species} column using their complete
#' scientific names. This is determined from lineage position rather than a
#' fixed list of rank labels: whenever the queried node is below a species
#' node, its own scientific name replaces the ancestral species name.
#'
#' If the database contains the \code{merged} table created by
#' \code{build_ncbi_taxonomy_database()}, obsolete tax IDs are automatically redirected to
#' their current NCBI tax IDs. The current ID is returned in
#' \code{query_taxid}; set \code{keep_input_taxid = TRUE} to retain the
#' original ID separately. If \code{name_aliases} is available, taxon-name
#' queries also recognize NCBI synonyms while returning the current tax ID and
#' current scientific name. If \code{species_equivalents} is available,
#' species names qualified with \code{s.l.} or \code{sensu lato} are resolved
#' to an unambiguous canonical binomial species. This makes sample-specific
#' records such as \code{Gnetum hainanense s.l. Won 511} return the same species
#' and tax ID as \code{Gnetum hainanense}, without collapsing named varieties,
#' subspecies, or subgroups.
#'
#' Bracketed provisional genera are displayed without brackets only when NCBI
#' supplies the corresponding unbracketed synonym for the same tax ID. For
#' example, tax ID 2995247 is returned as
#' \code{Polygonum chinense var. procumbens}; the \code{varietas} node is not
#' redirected to its parent species.
#'
#' A maximum recursion depth (`max_depth`) is used to avoid infinite loops
#' caused by malformed taxonomy entries.
#'
#' @return
#' A data frame containing:
#' \describe{
#'   \item{\code{query_taxid}}{
#'     Current NCBI tax ID after resolving obsolete IDs through
#'     \code{merged.dmp} and qualified species through
#'     \code{species_equivalents}.
#'   }
#'   \item{Taxonomic ranks}{
#'     One column per rank specified in \code{rank_name}.
#'   }
#'   \item{Matched name information}{
#'     When \code{tax_names} is used, additional columns reporting the original
#'     taxon names and the matched tax IDs are included.
#'   }
#' }
#' Missing ranks are filled with \code{NA}.
#'
#' @section Database Requirements:
#' The SQLite database must include tables structured identically to the
#' official NCBI taxonomy dump (`nodes.dmp`, `names.dmp`) after being imported
#' into SQLite. Required fields include parent tax IDs and rank annotations.
#' A \code{merged} table is optional for databases created by older package
#' versions, but is required to resolve obsolete tax IDs.
#'
#' @examples
#' \dontrun{
#' db <- "reference_databases/ncbi_taxonomy/taxonomy.sqlite"
#'
#' # Step 1: create the database once
#' build_ncbi_taxonomy_database(
#'   taxonomy = TRUE,
#'   accession = TRUE,
#'   taxdump_dir = "reference_databases/ncbi_taxonomy",
#'   sql_path = db
#' )
#'
#' # Step 2a: query taxonomy directly by taxid
#' lineage_by_id <- get_taxonomy(
#'   sql_path = db,
#'   tax_ids = c(9606, 3702)
#' )
#'
#' # Step 2b: query taxonomy directly by scientific name
#' lineage_by_name <- get_taxonomy(
#'   sql_path = db,
#'   tax_names = c("Homo sapiens", "Arabidopsis thaliana")
#' )
#'
#' # NCBI synonyms resolve to the current taxid and scientific name
#' current_name <- get_taxonomy(
#'   sql_path = db,
#'   tax_names = "Acer takesimense"
#' )
#'
#' # Step 2c: accession -> taxid -> taxonomy
#' accession_taxids <- get_taxid(
#'   sql_path = db,
#'   accessions = c("NC_045512.2", "NC_000001.11")
#' )
#' accession_lineages <- get_taxonomy(
#'   sql_path = db,
#'   tax_ids = accession_taxids$taxid
#' )
#' }
#'
#' @seealso
#' \code{\link{build_ncbi_taxonomy_database}} for creating or selectively updating the
#' required SQLite database, and \code{\link{get_taxid}} for converting
#' accession identifiers to tax IDs before lineage lookup.
#'
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbConnect dbExecute dbDisconnect dbWriteTable dbQuoteLiteral dbGetQuery
#' @importFrom dplyr left_join bind_cols filter rename join_by select  mutate arrange
#' @importFrom tidyr pivot_wider
#' @return query table
#' @export
get_taxonomy <- function(sql_path, tax_ids = NULL, tax_names = NULL, dereplicate = TRUE, rank_name = c("kingdom", "phylum","class","order","family", "genus" ,"species" ), show_elapsed = T, keep_input_taxid = FALSE) {
  start_time <- Sys.time()
  show_time("Starting taxonomic query...")
  on.exit(show_time("Done."), add = TRUE)
  if(show_elapsed) on.exit(show_time(secs = as.numeric(Sys.time() - start_time, units = "secs"), mode = 1), add = TRUE)
  # setting max_depth to search to avoid endless loop.
  max_depth = 100
  #
  con <- dbConnect(SQLite(),sql_path)
  # close connection when leaving function
  on.exit(dbDisconnect(con), add = TRUE)
  #
  skip_id_O <- F
  # Avoid inputting both parameters at the same time.
  stopifnot(xor(is.null(tax_ids), is.null(tax_names)))
  # raw order of ids

  if(!is.null(tax_ids)){
    id_raw_seq <- data.frame(
      input_taxid = tax_ids,
      query_taxid = .resolve_ncbi_taxids(con, tax_ids)
    )
    tax_ids <- unique(id_raw_seq$query_taxid)
  }
  # raw name and get their tax id
  if (!is.null(tax_names)) {
    skip_id_O <- T
    name_table <- if ("name_aliases" %in% DBI::dbListTables(con)) {
      "name_aliases"
    } else {
      "names"
    }
    name_raw <- DBI::dbGetQuery(
      con,
      paste0(
        "SELECT DISTINCT tax_id, name FROM ", name_table, " WHERE name IN (",
        paste(DBI::dbQuoteLiteral(con, tax_names), collapse = ","),
        ")"
      )
    )
    name_raw$tax_id <- .resolve_ncbi_taxids(con, name_raw$tax_id)
    name_raw <- distinct(name_raw, tax_id, name, .keep_all = TRUE)
    tax_ids <- name_raw$tax_id
    name_raw <- data.frame(query_tax_name = tax_names) %>% left_join(name_raw, by = join_by(query_tax_name == name))

    if (length(tax_ids) == 0) {
      n <- length(tax_names)

      tax_tbl <- as.data.frame(
        matrix(NA, nrow = n, ncol = length(rank_name),
               dimnames = list(NULL, rank_name))) %>%
        bind_cols(
          query_tax_name = tax_names,
          query_taxid = rep(NA_integer_, n),
          .
        )
      return(tax_tbl)
    }
  }

  # make a query ---------
  id_sql <- paste(DBI::dbQuoteLiteral(con, tax_ids), collapse = ",")


  q <- paste(
    "WITH RECURSIVE lineage(seed, depth, tax_id, parent_tax_id) AS (",
    "  SELECT tax_id AS seed, 0 AS depth, tax_id, parent_tax_id",
    "  FROM nodes WHERE tax_id IN (", id_sql, ")",
    "  UNION ALL",
    "  SELECT l.seed, l.depth + 1, p.tax_id, p.parent_tax_id",
    "  FROM lineage l",
    "  JOIN nodes p ON p.tax_id = l.parent_tax_id",
    "  WHERE p.parent_tax_id IS NOT NULL",
    "    AND p.tax_id <> p.parent_tax_id",
    "    AND l.depth < ", max_depth,
    ")",
    "SELECT l.seed AS query_taxid, l.depth,",
    "       n.tax_id, n.parent_tax_id, n.rank,",
    "       n.embl_code, n.division_id,",
    "       n.inherited_div_flag, n.genetic_code_id,",
    "       n.inherited_gc_flag, n.mito_genetic_code_id,",
    "       n.inherited_mgc_flag, n.genbank_hidden_flag,",
    "       n.subtree_root_flag, n.comments,",
    "       nm.name AS scientific_name",
    "FROM lineage l",
    "JOIN nodes n ON n.tax_id = l.tax_id",
    "LEFT JOIN names nm ON nm.tax_id = n.tax_id",
    "ORDER BY l.seed, l.depth;",
    sep = "\n"
  )

  query_res <- DBI::dbGetQuery(con, q) %>%
    select(query_taxid, depth, parent_tax_id, rank, scientific_name) %>%
    group_by(query_taxid) %>%
    mutate(
      rank = ifelse(
        depth == 0L & any(rank == "species" & depth > 0L),
        "species",
        rank
      )
    ) %>%
    ungroup() %>%
    filter(rank %in% rank_name) %>%
    arrange(query_taxid, depth) %>%
    distinct(query_taxid, rank, .keep_all = TRUE) %>%
    mutate(rank = rank %>% factor(levels = rev(rank_name))) %>%
    arrange(rank) %>%
    pivot_wider(names_from = rank, id_cols = query_taxid, values_from =  scientific_name, values_fill = NA)

  missing_ranks <- setdiff(rank_name, names(query_res))
  query_res[missing_ranks] <- NA_character_
  query_res <- query_res[, c("query_taxid", rank_name), drop = FALSE]

  if(!is.null(tax_names)){
    final_query_res <- left_join(name_raw, query_res, by = join_by(tax_id == query_taxid)) %>%
      dplyr::rename(query_taxid = tax_id)
  }
  if(!skip_id_O){
    final_query_res <- left_join(
      id_raw_seq,
      query_res,
      by = join_by(query_taxid)
    )
    if (!isTRUE(keep_input_taxid)) {
      final_query_res <- final_query_res %>% select(-input_taxid)
    }
  }
  # remove replicated results -------------------------
  if(dereplicate){
    final_query_res <- final_query_res %>% distinct(query_taxid, .keep_all = T)
  }
    return(final_query_res)

  show_time("Done.")

}

# Query taxid by accesion number ------------------------------
#' Retrieve NCBI taxonomic identifiers (taxid) for accessions
#'
#' Query an SQLite accession-to-taxid database to obtain taxonomic identifiers
#' corresponding to a vector of NCBI accessions. Both versioned accessions
#' (e.g. \code{AB123456.1}) and unversioned accessions (e.g. \code{AB123456})
#' are supported and are detected automatically.
#'
#' The function returns results in the same order as the input accessions.
#' Accessions not found in the database will have \code{NA} as their taxid.
#'
#' @param sql_path Character scalar. Path to the SQLite database generated by
#'   \code{build_ncbi_taxonomy_database()}, containing the \code{accession2taxid} table.
#' @param accessions Character vector of NCBI accessions to be queried.
#'   Versioned identifiers (for example, \code{"AB123456.1"}) and unversioned
#'   identifiers (for example, \code{"AB123456"}) may be mixed in one call.
#' @param show_elapsed Logical. Whether to display elapsed time information
#'   using \code{show_time()}. Default is \code{TRUE}.
#'
#' @return A \code{data.frame} with two columns:
#' \describe{
#'   \item{accession}{Input accession identifiers, in the same order as provided.}
#'   \item{taxid}{NCBI taxonomic identifier corresponding to each accession.}
#' }
#'
#' @details
#' The function checks whether the database contains a table named
#' \code{accession2taxid}. If not found, execution is stopped with an error.
#' Inputs ending in a dot followed by digits are queried against
#' \code{accession.version}; all other inputs are queried against
#' \code{accession}. Results retain the original input values and order.
#' When the database contains a \code{merged} table, obsolete accession tax IDs
#' are automatically replaced by their current NCBI tax IDs before return.
#' SQL queries are constructed internally using \code{IN (...)} clauses;
#' therefore, extremely large accession vectors may be limited by SQLite
#' query length constraints.
#'
#' @examples
#' \dontrun{
#' db <- "reference_databases/ncbi_taxonomy/taxonomy.sqlite"
#'
#' # Step 1: build taxonomy and accession tables
#' build_ncbi_taxonomy_database(
#'   taxonomy = TRUE,
#'   accession = TRUE,
#'   taxdump_dir = "reference_databases/ncbi_taxonomy",
#'   sql_path = db
#' )
#'
#' # Step 2: map mixed versioned and unversioned accessions to taxids
#' accession_taxids <- get_taxid(
#'   sql_path = db,
#'   accessions = c("NC_045512.2", "NC_000001")
#' )
#'
#' # Step 3: retrieve the corresponding taxonomic lineages
#' accession_lineages <- get_taxonomy(
#'   sql_path = db,
#'   tax_ids = accession_taxids$taxid
#' )
#' }
#'
#' @seealso
#' \code{\link{build_ncbi_taxonomy_database}} for creating or selectively updating the
#' accession-to-taxid database, and \code{\link{get_taxonomy}} for retrieving
#' the lineage associated with returned tax IDs.
#'
#' @import DBI
#' @import RSQLite
#' @importFrom dplyr left_join rename
#' @importFrom stringr str_detect
#'
#' @export
get_taxid <- function(sql_path, accessions, show_elapsed = T){
  show_time("Starting taxid query...")
  start_time <- Sys.time()
  con <- dbConnect(SQLite(),sql_path)

  if(show_elapsed) on.exit(show_time(secs = as.numeric(Sys.time() - start_time, units = "secs"), mode = 1))
  on.exit(dbDisconnect(con), add = T)
  # test <- tbl(con,  "accession2taxid")  %>% select(accession.version, taxid) %>%  filter(accession.version %in% accessions) %>% show_query()

  if(!"accession2taxid" %in% dbListTables(con)){
    stop(
      "The accession2taxid table is missing. Build it with build_ncbi_taxonomy_database().",
      call. = FALSE
    )
  }
  if(!is.character(accessions)){
    stop("accessions must be a character vector.", call. = FALSE)
  }
  if(length(accessions) == 0L){
    return(data.frame(accession = character(), taxid = integer()))
  }

  query_accessions <- function(values, column){
    values <- unique(values[!is.na(values)])
    if(length(values) == 0L){
      return(data.frame(accession = character(), taxid = integer()))
    }
    quoted_values <- paste(DBI::dbQuoteLiteral(con, values), collapse = ",")
    query <- sprintf(
      paste(
        "SELECT `%s` AS accession, taxid",
        "FROM accession2taxid",
        "WHERE `%s` IN (%s)"
      ),
      column, column, quoted_values
    )
    DBI::dbGetQuery(con, query)
  }

  is_versioned <- !is.na(accessions) & grepl("\\.[0-9]+$", accessions)
  versioned_res <- query_accessions(
    accessions[is_versioned],
    "accession.version"
  )
  unversioned_res <- query_accessions(
    accessions[!is_versioned & !is.na(accessions)],
    "accession"
  )
  res <- rbind(versioned_res, unversioned_res)
  if (nrow(res) > 0L) {
    res$taxid <- .resolve_ncbi_taxids(con, res$taxid)
  }

  show_time("Done!")
  data.frame(
    accession = accessions,
    taxid = res$taxid[match(accessions, res$accession)]
  )
}


# Unzip the blast database --------------------------------------
#' Decompress Downloaded NCBI BLAST Database Fragments
#'
#' This function extracts `.tar.gz` BLAST database fragments previously
#' downloaded from the NCBI BLAST FTP server. Each compressed archive is
#' decompressed into a specified directory, producing the files required for
#' local BLAST database use. Parallel processing is supported via the
#' \pkg{future} and \pkg{foreach} frameworks.
#'
#' @param file_path (character) Directory containing downloaded BLAST
#'   `.tar.gz` fragment files. Must be an absolute path to avoid
#'   extraction errors. Default is the current working directory.
#'
#' @param database_path (character) Directory where decompressed BLAST
#'   database files will be written. The directory must be writable.
#'
#' @param threads Integer. Number of parallel workers used for decompression.
#'   If `NULL`, decompression is performed sequentially. Default is `NULL`.
#'
#' @param show_elapsed Logical. If `TRUE`, the total elapsed time is reported
#'   using the internal \code{show_time()} function. Default is `TRUE`.
#'
#' @details
#' The function scans \code{file_path} for all files matching the pattern
#' \code{"tar.gz"} and extracts each archive using \code{untar()}.
#' When \code{threads} is specified, parallel processing is performed using
#' \pkg{future} and \pkg{foreach} (with \code{\%dofuture\%}).
#'
#' Each `.tar.gz` fragment typically corresponds to a portion of a BLAST
#' database (e.g., `nt.00.tar.gz`, `nt.01.tar.gz`). After extraction, the
#' resulting files (e.g., `.nhr`, `.nin`, `.nsq`) form a complete BLAST
#' database when all fragments are included.
#'
#' Future plans are automatically restored to sequential mode upon exit.
#'
#' @return
#' Invisibly returns \code{NULL}. Decompressed BLAST database files are
#' written into \code{database_path}.
#'
#' @section Notes:
#' \itemize{
#'   \item All paths should be absolute to avoid expansion issues during
#'         extraction.
#'   \item Large databases such as \code{nt} may require substantial disk
#'         space after decompression.
#'   \item This function does not perform MD5 checksum verification.
#' }
#'
#' @examples
#' \dontrun{
#' # Decompress all downloaded BLAST fragments
#' prepare_blast_database(
#'   file_path = "reference_databases/blast_downloads",
#'   database_path = "reference_databases/blast_database",
#'   threads = 8
#' )
#' }
#'
#' @import foreach
#' @importFrom utils untar
#' @importFrom future plan multisession sequential
#' @importFrom doFuture `%dofuture%`
#' @export
prepare_blast_database <- function(file_path = ".", database_path = ".", threads = NULL, show_elapsed = T){
  start_time <- Sys.time()
  if(show_elapsed) on.exit(show_time(secs = as.numeric(Sys.time() - start_time, units = "secs"), mode = 1), add = TRUE)
  show_time("Starting to decompress the files...")
  show_time("Unzipping data...")
  on.exit(plan(sequential), add = T)
  if (is.null(threads)) {
    plan(sequential)
  } else {
    plan(multisession, workers = threads)
  }

  file_l <- list.files(file_path, full.names = TRUE, pattern = "tar.gz$")

  # handlers(global = TRUE)
  # handlers("progress")

  #with_progress({

   # p <- progressor(along = 1:length(file_l))

    foreach(fl = seq_along(file_l)) %dofuture% {

      untar(
        tarfile = file_l[fl],
        exdir = database_path
      )

    #  p()
    }
  #})


  show_time("Done.")
}



## LCA function ---------------------------------
#' Infer Last Common Ancestor (LCA) taxonomy from BLAST tabular results
#'
#' This function infers a consensus taxonomic assignment for each query
#' sequence based on BLAST tabular output. For each query, only the best hits
#' (defined by maximum identity, maximum query coverage, and minimum e-value)
#' are retained. Taxonomic ranks are then evaluated independently, and a rank
#' is assigned only if all retained hits share the same value at that level.
#' The deepest rank with a non-ambiguous consensus is reported as the inferred
#' LCA rank.
#'
#' @param df (data.frame or tibble) BLAST tabular results containing query
#'   identifiers, alignment statistics, and taxonomic annotations. The table
#'   must include at least columns for query ID, percent identity, query
#'   coverage, e-value, taxonomy ranks, taxonomic identifiers, and subject
#'   accessions.
#'
#' @param seqs_ID (character) column name in \code{df} identifying query
#'   sequences (e.g., \code{"qseqid"}).
#'
#' @param tax_levels (character vector) ordered taxonomic ranks used for LCA
#'   inference, from higher to lower rank (e.g.,
#'   \code{c("kingdom","phylum","class",
#'   "order","family","genus","species")}).
#'
#' @param taxid_colum (character) column name containing NCBI taxonomic
#'   identifiers associated with BLAST hits (e.g., \code{"staxids"}).
#'
#' @param accession_colum (character) column name containing subject accession
#'   identifiers from BLAST hits (e.g., \code{"sseqid"}).
#'
#' @param identity_colum (character) column name for BLAST percent identity
#'   values (e.g., \code{"pident"}).
#'
#' @param coverage_colum (character) column name for BLAST query coverage
#'   values (e.g., \code{"qcovs"}).
#'
#' @details
#' The function proceeds as follows:
#' \itemize{
#'   \item Groups BLAST hits by query sequence.
#'   \item Retains only hits with:
#'     \itemize{
#'       \item maximum percent identity,
#'       \item maximum query coverage,
#'       \item minimum e-value.
#'     }
#'   \item For each taxonomic rank, assigns a value only if all retained hits
#'         share the same annotation at that rank; otherwise assigns \code{NA}.
#'   \item Determines the deepest (most specific) rank with a non-\code{NA}
#'         consensus as the inferred LCA rank.
#'   \item Concatenates taxonomic IDs and accession IDs across retained hits
#'         for traceability.
#'   \item Generates a \code{scientific_name} column, appending \code{"sp."}
#'         identifiers when the consensus rank is above species level.
#' }
#'
#' This approach corresponds to a *rank-wise consensus LCA* based on best BLAST
#' hits, rather than a phylogenetic tree-based LCA.
#'
#' @return
#' A tibble with one row per query sequence, including:
#' \itemize{
#'   \item consensus taxonomic annotations at each rank,
#'   \item inferred LCA rank (\code{rank}),
#'   \item constructed scientific name (\code{scientific_name}),
#'   \item concatenated taxonomic IDs (\code{taxid}),
#'   \item concatenated accession identifiers (\code{accession}),
#'   \item representative identity and coverage values.
#' }
#'
#' @importFrom dplyr coalesce sym syms summarise filter group_by select mutate
#' @importFrom stats na.omit
#'
#' @examples
#' blast_hits <- data.frame(
#'   qseqid = c("query_1", "query_1", "query_2"),
#'   sseqid = c("AB000001.1", "AB000002.1", "AB000003.1"),
#'   staxids = c(9606, 63221, 10090),
#'   pident = c(99.5, 99.5, 100),
#'   qcovs = c(100, 100, 100),
#'   evalue = c(1e-40, 1e-40, 1e-50),
#'   kingdom = c("Eukaryota", "Eukaryota", "Eukaryota"),
#'   phylum = c("Chordata", "Chordata", "Chordata"),
#'   class = c("Mammalia", "Mammalia", "Mammalia"),
#'   order = c("Primates", "Primates", "Rodentia"),
#'   family = c("Hominidae", "Hominidae", "Muridae"),
#'   genus = c("Homo", "Homo", "Mus"),
#'   species = c("Homo sapiens", "Homo sp.", "Mus musculus")
#' )
#'
#' lca_table <- find_lca(
#'   blast_hits,
#'   tax_levels = c("kingdom", "phylum", "class", "order", "family",
#'                  "genus", "species")
#' )
#' @section Citations: Wang, Y., Korneliussen, T. S., Holman, L. E., Manica, A., & Pedersen, M. W. (2022). ngsLCA - A toolkit for fast and flexible lowest common ancestor inference and taxonomic profiling of metagenomic data. Methods in Ecology and Evolution, 13, 2699-2708. \url{https://doi.org/10.1111/2041-210X.14006}
#' @export
find_lca <- function(df, seqs_ID = "qseqid", tax_levels = c("kingdom", "phylum","class",
                                                            "order","family", "genus" ,"species" ),
                     taxid_colum = "staxids", accession_colum = "sseqid",
                     identity_colum = "pident", coverage_colum = "qcovs" ) {

  get_consensus <- function(x) {
    vals <- unique(na.omit(x))
    if (length(vals) == 1) vals else NA_character_
  }

  lca_raw <- df %>%
    group_by(!!sym(seqs_ID)) %>%
    filter(pident ==max(pident), evalue == min(evalue), qcovs == max(qcovs )) %>%
    summarise(
      across(all_of(tax_levels), get_consensus),
      taxid = paste(!!sym(taxid_colum), collapse = ";"),
      accession = paste(!!sym(accession_colum), collapse = ";"),
      coverage = (!!sym(coverage_colum) %>% unique)[1],
      identity = (!!sym(identity_colum) %>% unique)[1],
      .groups = "drop"
    )

    t_rank <- lca_raw %>% select(all_of(rev(tax_levels))) %>% apply(1, function(x,tax_levels){rev(tax_levels)[which(!(x %>% is.na))[1]]},tax_levels = tax_levels)
    summ_df <- lca_raw %>%
      mutate(rank = t_rank) %>%
      group_by(across(all_of(tax_levels))) %>%
      mutate(scientific_name = coalesce(!!! syms(rev(tax_levels))),
             scientific_name = case_when(
               rank != tax_levels[length(tax_levels)] ~ paste(scientific_name, "sp.", dplyr::row_number()),
               .default = scientific_name
             ) )

  summ_df %>% return()

}
