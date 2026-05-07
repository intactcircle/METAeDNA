# install.packages(c("tidyverse", "DBI", "RSQLite"))
# library(tidyverse)
# library(DBI)
# library(RSQLite)
# library(curl)
# Constructing the SQL database for taxdump and accession2id.-------------------------------------------------------------
#' Construct an SQLite database for NCBI taxdump and accession-to-taxid mappings
#'
#' This function builds an SQLite database from NCBI taxonomy resources. It can
#' optionally download and parse the NCBI taxdump archive (nodes and names) and
#' the large \code{accession2taxid} mapping files, then write them into an
#' SQLite database for fast querying by downstream functions such as
#' \code{get_taxid()} and \code{get_taxonomy()}. The function supports
#' incremental construction of taxonomy-only, accession-only, or combined
#' databases, and can overwrite an existing database if requested.
#'
#' @param taxonomy (logical) if \code{TRUE}, create taxonomy tables from
#'   \code{nodes.dmp} and \code{names.dmp} in the NCBI taxdump archive. These
#'   tables are written as \code{nodes} and \code{names} in the SQLite database.
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
#' @param taxdump_dir (character) directory used to store downloaded taxdump
#'   and accession-to-taxid files. This directory must have sufficient free
#'   disk space (typically at least 10–60 GB, depending on whether accession
#'   tables and WGS data are included).
#'
#' @param sql_path (character) file path where the SQLite database will be
#'   created. Using an absolute path is recommended. The target file system
#'   should provide at least 60 GB of free space when building a full
#'   accession2taxid database with WGS data.
#'
#' @param overwrite (logical) if \code{TRUE}, remove any existing database file
#'   at \code{sql_path} before creating a new one. If \code{FALSE} and the file
#'   exists, the function emits a warning and does not overwrite the database.
#'
#' @param download (logical) if \code{TRUE}, download the required taxdump and
#'   accession2taxid files from the NCBI FTP site into \code{taxdump_dir}. If
#'   \code{FALSE}, the function assumes that the relevant files already exist
#'   in \code{taxdump_dir} with their original NCBI file names.
#'
#' @param resume (logical) if \code{TRUE}, attempt to resume interrupted
#'   downloads using \code{curl::multi_download()}.
#'
#' @param create_index (logical) if \code{TRUE}, create indexes on key columns
#'   (e.g., \code{parent_tax_id} in \code{nodes}, \code{tax_id} in \code{names},
#'   and \code{accession} or \code{taxid} in \code{accession2taxid}) to improve
#'   query performance at the cost of additional disk space and initial build
#'   time.
#'
#' @details
#'
#' Building the accession-to-taxid table, especially with WGS data, requires
#' substantial disk space and I/O time. Progress messages and approximate
#' record counts are printed during import.
#'
#' @return
#' Invisibly returns \code{NULL}. The main result is the SQLite database file
#' created at \code{sql_path}, containing some or all of the following tables:
#' \code{nodes}, \code{names}, and \code{accession2taxid}, depending on the
#' arguments supplied.
#' @importFrom DBI dbConnect dbExecute dbDisconnect dbWriteTable
#' @importFrom RSQLite SQLite
#' @importFrom readr read_delim read_tsv
#' @importFrom curl curl_download
#' @importFrom utils untar
#' @importFrom stringr str_remove
#'
#' @examples
#' \dontrun{
#' # Build taxonomy-only database in the current directory
#' make_database(
#'   taxonomy   = TRUE,
#'   accession  = FALSE,
#'   include_WGS = FALSE,
#'   taxdump_dir = ".",
#'   sql_path    = "./taxdump.sqlite",
#'   overwrite   = FALSE
#' )
#'
#' # Build taxonomy + accession2taxid (GenBank) database with overwrite
#' make_database(
#'   taxonomy   = TRUE,
#'   accession  = TRUE,
#'   include_WGS = FALSE,
#'   taxdump_dir = "/data/ncbi_taxdump",
#'   sql_path    = "/data/ncbi_taxdump/taxdump.sqlite",
#'   overwrite   = TRUE
#' )
#' }
#' @export
make_database <- function(taxonomy = TRUE, accession = FALSE, include_WGS = FALSE, taxdump_dir = ".", sql_path = "./taxdump.sqlite", overwrite = F, download = T,resume = T,create_index = T){
  start_time <- Sys.time()
  on.exit(show_time(sec = as.numeric(Sys.time() - start_time , units = "secs") , mode =1))
  #browser()
  message("[Warning]: Please ensure that you use the complete absolute path.")
  show_time("Start constructing the SQL database...")
  if(file.exists(sql_path)&!overwrite){
    message("[Warning]: The database already exists. You can set overwrite = TRUE to replace it.")
  }
  if(file.exists(sql_path)&overwrite){
    unlink(sql_path)
  }
  if(taxonomy){
    #
    if(download){
      show_time("Download taxdump file...")
      curl::multi_download("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", file.path(taxdump_dir, "taxdump.tar.gz"), resume = resume)
    }

    untar( file.path(taxdump_dir, "taxdump.tar.gz"), exdir = taxdump_dir)

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
    names <- read_delim(file.path(taxdump_dir, "names.dmp"), delim = "\t|\t" ,col_names = names_cols, show_col_types = F, progress = F) %>%
      mutate(name_class = name_class %>% str_remove("\t[|]$")) %>%
      filter(name_class == "scientific name") %>%
      select(tax_id = tax_id, name = name_txt)

    # 3) write SQL table -------------------------------------------
    con <- dbConnect(SQLite(), sql_path)
    # close connection when leaving function

    #
    dbExecute(con, "PRAGMA journal_mode=DELETE;")
    # remove old data if existed.
    dbExecute(con, "DROP TABLE IF EXISTS nodes;")
    dbExecute(con, "DROP TABLE IF EXISTS names;")
    dbExecute(con, "DROP TABLE IF EXISTS merged;")
    # write table
    dbWriteTable(con, "nodes", nodes, overwrite = TRUE)
    dbWriteTable(con, "names", names, overwrite = TRUE)
    #
    show_time("Indexing [taxonomy] table ...")
    # 4) indexing ----------------
    dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_nodes_parent ON nodes(parent_tax_id);")
    dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_names_tax ON names(tax_id);")
    dbDisconnect(con)
  }
  # ----------------------
  if(include_WGS & download){
    show_time("Download [nucl_wgs.accession2taxid.gz] and [nucl_gb.accession2taxid.gz] files...")
    curl::multi_download(c("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
                           "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"),
                         destfile = file.path(taxdump_dir,c("nucl_wgs.accession2taxid.gz","nucl_gb.accession2taxid.gz")),
                         resume = resume)
  }else if(download){
    show_time("Download [nucl_gb.accession2taxid.gz] file...")
    curl::multi_download("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
                         destfile = file.path(taxdump_dir, "nucl_gb.accession2taxid.gz"), resume = resume)
  }
  # ---------------------------
  wr_print <- T
  if(accession){
    show_time("Constructing the accession2id table...")
    message("[Warning]: Since you have chosen to build the accession-to-taxid database, please ensure that you have more than 60 GB of available disk space, as this database is extremely large.")
    message("[Warning]: You may also download the file nucl_gb.accession2taxid.gz manually and place it in the current directory, but please do not modify its filename.")
    message("[INFO]: You may enjoy a 20-minute coffee break; the exact duration depends on your network speed. -.-!")
    wr_print <- F


    #
    con <- dbConnect(SQLite(), sql_path)
    #
    dbExecute(con, "DROP TABLE IF EXISTS accession2taxid;")
    # 1) file connection ----
    gz_con <- gzfile( file.path(taxdump_dir, "nucl_gb.accession2taxid.gz"), open = "rb" )
    # header
    header_line <- readLines(gz_con,n = 1)
    header <- read_tsv(I(header_line), show_col_types = F) %>% colnames
    show_time(sprintf("Processing accession2taxid records..."), mode = 0)
    # 2) chunked -------------------------------------
    n_proccessed_rec <- 0
    repeat{
      chunked_data <- readLines(gz_con, n = 1000000)
      if(length(chunked_data) == 0 ){break}
      n_proccessed_rec <- n_proccessed_rec + length(chunked_data)
      dbWriteTable(con, "accession2taxid" ,read_tsv(I(chunked_data), col_names = header,show_col_types=F,  progress = F), append = TRUE )
      show_time(sprintf("Processing records in [nucl_gb.accession2taxid.gz]: %d...\r", n_proccessed_rec), mode = 2)
    }
    cat("\n")

    close(gz_con)
    dbDisconnect(con)
  }
  # ------------------------------
  if(include_WGS){
    con <- dbConnect(SQLite(), sql_path)
    if(wr_print){
      show_time("Constructing the accession2id table...")
      message("[Warning]: Since you have chosen to build the accession-to-taxid database, please ensure that you have more than 60 GB of available disk space, as this database is extremely large.")
      message("[Warning]: You may also download the file nucl_gb.accession2taxid.gz manually and place it in the current directory, but please do not modify its filename.")
      message("[INFO]: You may enjoy a 20-minute coffee break; the exact duration depends on your network speed. -.-!")

    }
        # 1) file connection ----
      gz_con <- gzfile( file.path(taxdump_dir, "nucl_wgs.accession2taxid.gz"), open = "rb" )
      # header
      header_line <- readLines(gz_con,n = 1)
      header <- read_tsv(I(header_line), show_col_types = F) %>% colnames
      show_time(sprintf("Processing WGS records..."), mode = 0)
      # 2) chunked -------------------------------------
      n_proccessed_rec <- 0
      repeat{
        chunked_data <- readLines(gz_con, n = 1000000)
        if(length(chunked_data) == 0 ){break}
        n_proccessed_rec <- n_proccessed_rec + length(chunked_data)
        dbWriteTable(con, "accession2taxid" ,read_tsv(I(chunked_data), col_names = header,show_col_types= F,  progress = F), append = TRUE )
        show_time(sprintf("Processing records in [nucl_wgs.accession2taxid.gz]: %d...\r", n_proccessed_rec), mode = 2)
      }
      cat("\n")
      close(gz_con)
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
#'
#' @details
#' The function constructs a recursive SQL query using a Common Table
#' Expression (CTE) to follow taxonomic lineages upward from each query tax ID.
#' Returned results include only user-defined ranks. When `tax_names` is used,
#' the function also returns the matched tax IDs for each supplied name.
#'
#' A maximum recursion depth (`max_depth`) is used to avoid infinite loops
#' caused by malformed taxonomy entries.
#'
#' @return
#' A data frame containing:
#' \describe{
#'   \item{\code{query_tax_id}}{
#'     Input tax ID or the tax ID matched to the given scientific name.
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
#'
#' @examples
#' \dontrun{
#' # Example: query by tax IDs
#' get_taxonomy(
#'   sql_path = "ncbi_taxonomy.sqlite",
#'   tax_ids = c(9606, 3702)
#' )
#'
#' # Example: query by scientific names
#' get_taxonomy(
#'   sql_path = "ncbi_taxonomy.sqlite",
#'   tax_names = c("Homo sapiens", "Arabidopsis thaliana")
#' )
#' }
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbConnect dbExecute dbDisconnect dbWriteTable dbQuoteLiteral dbGetQuery
#' @importFrom dplyr left_join bind_cols filter rename join_by select  mutate arrange
#' @importFrom tidyr pivot_wider
#' @return query table
#' @export
get_taxonomy <- function(sql_path, tax_ids = NULL, tax_names = NULL, dereplicate = TRUE, rank_name = c("kingdom", "phylum","class","order","family", "genus" ,"species" ), show_elapsed = T) {
  start_time <- Sys.time()
  show_time("Starting taxonomic query...")
  on.exit( show_time("Done."))
  if(show_elapsed) on.exit(show_time(sec = as.numeric(Sys.time() - start_time , units = "secs") , mode =1), add = T)
  # setting max_depth to search to avoid endless loop.
  max_depth = 100
  #
  # browser()
  con <- dbConnect(SQLite(),sql_path)
  # close connection when leaving function
  on.exit(dbDisconnect(con))
  #
  skip_id_O <- F
  # Avoid inputting both parameters at the same time.
  stopifnot(xor(is.null(tax_ids), is.null(tax_names)))
  # raw order of ids

  if(!is.null(tax_ids)){
    id_raw_seq <- data.frame(query_tax_id = tax_ids)
  }
  # raw name and get their tax id
  if (!is.null(tax_names)) {
    skip_id_O <- T
    name_raw <- DBI::dbGetQuery(
      con,
      paste0(
        "SELECT DISTINCT tax_id, name FROM names WHERE name IN (",
        paste(DBI::dbQuoteLiteral(con, tax_names), collapse = ","),
        ")"
      )
    )
    tax_ids <- name_raw$tax_id
    name_raw <- data.frame(query_tax_name = tax_names) %>% left_join(name_raw, by = join_by(query_tax_name == name))

    if (length(tax_ids) == 0) {
      n <- length(tax_names)

      tax_tbl <- as.data.frame(
        matrix(NA, nrow = n, ncol = length(rank_name),
               dimnames = list(NULL, rank_name))) %>%
        bind_cols(query_tax_name = tax_names,.)
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
    "SELECT l.seed AS query_tax_id, l.depth,",
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
    select(query_tax_id, parent_tax_id, rank, scientific_name) %>%
    filter(rank %in% rank_name) %>%
    mutate(rank = rank %>% factor(levels = rev(rank_name))) %>%
    arrange(rank) %>%
    pivot_wider(names_from = rank, id_cols = query_tax_id, values_from =  scientific_name, values_fill = NA)
  #


  if(!is.null(tax_names)){
    final_query_res <- left_join(name_raw, query_res, by = join_by(tax_id == query_tax_id)) %>%
      dplyr::rename(query_tax_id = tax_id)
  }
  if(!skip_id_O){
    final_query_res <- left_join(id_raw_seq, query_res, by = join_by(query_tax_id))
  }
  # remove replicated results -------------------------
  if(dereplicate){
    final_query_res <- final_query_res %>% distinct(query_tax_id, .keep_all = T)
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
#'   \code{make_database()}, containing the \code{accession2taxid} table.
#' @param accessions Character vector of NCBI accessions to be queried.
#'   All accessions should be either versioned or unversioned.
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
#' SQL queries are constructed internally using \code{IN (...)} clauses;
#' therefore, extremely large accession vectors may be limited by SQLite
#' query length constraints.
#'
#' @examples
#' \dontrun{
#' db <- "accession2taxid.sqlite"
#' acc <- c("AB123456.1", "NM_001200025.1")
#' res <- get_taxid(db, acc)
#' head(res)
#' }
#'
#' @seealso \code{\link{make_database}}
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

  if(show_elapsed) on.exit(show_time(sec = as.numeric(Sys.time() - start_time , units = "secs") , mode =1))
  on.exit(dbDisconnect(con), add = T)
  # test <- tbl(con,  "accession2taxid")  %>% select(accession.version, taxid) %>%  filter(accession.version %in% accessions) %>% show_query()

  if(all(!dbListTables(con) == "accession2taxid")){
    stop("Database not found. Please use make_database() to build the accession-to-taxid database first.")
  }
  if(str_detect(accessions[1], "[.]")){
    q <- sprintf("SELECT  `accession.version`,  `taxid`
          FROM `accession2taxid`
          WHERE (`accession.version` IN (%s))",  paste(sprintf("'%s'", accessions), collapse = ",") )

    res <- DBI::dbGetQuery(con, q) %>% dplyr::rename(accession = `accession.version`)

      }else{
    q <- sprintf("SELECT  `accession`,  `taxid`
          FROM `accession2taxid`
          WHERE (`accession` IN (%s))",  paste(sprintf("'%s'", accessions), collapse = ",") )

    res <- DBI::dbGetQuery(con, q)

      }
  show_time("Done！")
  #

  data.frame(accession = accessions) %>%
    left_join(res,by = join_by(accession)) %>%
    return()


}


# Download BLAST database files ---------------------------------------------
#' Download selected NCBI BLAST databases
#'
#' This function downloads one or more predefined NCBI BLAST databases
#' (e.g., `nt`, `SSU_euk`, `ITS_fungal`) from the official NCBI FTP server.
#' The function retrieves the list of available `.tar.gz` BLAST database
#' fragments, identifies matching files for the selected database type,
#' and downloads them to a user-specified directory.
#'
#' @param db (character) One or more BLAST database names to download.
#'   Supported values include:
#'   \itemize{
#'     \item{\code{"nt"}}{ Nucleotide collection}
#'     \item{\code{"nt_euk"}}{ Eukaryotic subset of nt}
#'     \item{\code{"nt_vir"}}{ Viral subset of nt}
#'     \item{\code{"SSU_pro"}}{ Prokaryotic SSU rRNA}
#'     \item{\code{"SSU_euk"}}{ Eukaryotic SSU rRNA}
#'     \item{\code{"LSU_pro"}}{ Prokaryotic LSU rRNA}
#'     \item{\code{"LSU_euk"}}{ Eukaryotic LSU rRNA}
#'     \item{\code{"ITS_euk"}}{ Eukaryotic ITS}
#'     \item{\code{"SSU_fungal"}}{ Fungal SSU (18S)}
#'     \item{\code{"LSU_fungal"}}{ Fungal LSU (28S)}
#'     \item{\code{"ITS_fungal"}}{ Fungal ITS RefSeq}
#'   }
#'
#' @param database_path (character) Path to the directory where downloaded files
#'   will be saved. The directory will be created if it does not exist.
#'
#' @param show_elapsed (logical) If `TRUE`, the elapsed download time is
#'   displayed using `show_time()`. Default is `TRUE`.
#'
#' @details
#' The function identifies available BLAST database files by parsing the FTP
#' directory listing from:
#'
#' \verb{https://ftp.ncbi.nlm.nih.gov/blast/db/}
#'
#' Matching `.tar.gz` fragments are selected using regular expression filters.
#' All matched file names, metadata, and estimated sizes are written to a CSV
#' file (`database_file_info.csv`) inside `database_path`. The download is
#' performed using `multi_download()`, which supports resuming interrupted
#' downloads.
#'
#' File size values in the FTP listing are parsed as GB for estimation
#' purposes only and do not represent exact sizes.
#'
#' @return
#' Invisibly returns the result of `multi_download()`, typically a vector
#' containing destination file paths.
#'
#' @section Notes:
#' \itemize{
#'   \item The function downloads only original `.tar.gz` fragments; corresponding `.md5` checksum files are not downloaded.
#'   \item Only BLAST database files are retrieved; decompression and indexing must be performed separately using \code{blastdbcmd} or \code{update_blastdb.pl}.
#'   \item Large databases such as \code{nt} may require substantial storage and bandwidth.
#' }
#'
#' @examples
#' \dontrun{
#' # Download fungal ITS BLAST database
#' download_blast_db("ITS_fungal", database_path = "blastdb")
#'
#' # Download nt and nt_euk
#' download_blast_db(c("nt", "nt_euk"), database_path = "/data/blast/")
#' }
#' @importFrom dplyr filter mutate case_when select
#' @importFrom readr read_csv write_csv
#' @importFrom stringr str_detect str_split str_remove str_squish str_trim
#' @importFrom curl multi_download
#' @export
download_blast_db <- function( db = c("nt","nt_euk","nt_vir","SSU_pro","SSU_euk",
                                      "LSU_pro","LSU_euk","ITS_euk","SSU_fungal",
                                      "LSU_fungal", "ITS_fungal"),
                               database_path = ".", show_elapsed = T ) {
  if (!dir.exists(database_path)) dir.create(database_path, recursive = TRUE)
  start_time <- Sys.time()
  if(show_elapsed) on.exit(show_time(sec = Sys.time() - start_time , mode =1))
  show_time(paste0("Download NCBI [", db,"] BLAST database..."))
  base_url <- "https://ftp.ncbi.nlm.nih.gov/blast/db/"
  # Files list ----------------
  html <- tryCatch(readLines(base_url, warn = FALSE), error = function(e) "")
  file_info_table <- data.frame(
    file_name = html %>% str_remove("[<]a.href=\"") %>% str_remove("\".*"),
    meta_data = html %>% str_remove(".*/a[>]") %>% str_trim() %>% str_squish()
  ) %>%
    mutate(
      date = str_split(meta_data, " ") %>% sapply( `[`, 1),
      size = str_split(meta_data, " ") %>% sapply( `[`, 3),
    ) %>%
    filter(str_detect(file_name, "tar.gz")) %>%
    filter(!str_detect(file_name, "md5$")) %>%
    mutate(
      size_type = case_when(
        str_detect(size, "G") ~ "G",
        str_detect(size, "M") ~ "M",
        str_detect(size, "K") ~ "K",
        .default = "B",
      ),
      size = str_remove(size, "[A-Za-z].*") %>% as.numeric(),
      size = case_when(
        size_type == "G" ~ size,
        size_type == "M"~ size/1024,
        size_type == "K"~ size/(1024*1024),
        .default = size/(1024*1024*1024)
      ),
      url = paste0(base_url,file_name)
    ) %>%
    select(-size_type,-meta_data)
  # ---------------------------------
  download_type <- c(nt = "^nt[.]", nt_euk = "nt_euk", nt_vir = "nt_viruses",
                     SSU_pro = "16S_ribosomal_RNA", SSU_euk = "SSU_eukaryote_rRNA",
                     LSU_pro = "LSU_prokaryote_rRNA", LSU_euk = "LSU_eukaryote_rRNA",
                     ITS_euk = "ITS_eukaryote_sequences",
                     SSU_fungal = "18S_fungal_sequences",
                     LSU_fungal =  "28S_fungal_sequences",
                     ITS_fungal = "ITS_RefSeq_Fungi")

  target_table <- file_info_table %>% filter(str_detect(file_name,download_type[db]))
  target_table %>% write_csv(file.path(database_path,"database_file_info.csv"))

  target_table  <-  read_csv(file.path(database_path,"database_file_info.csv"), show_col_types = F)
  message("[INFO] Estimated file size: " %>% paste0(target_table$size %>% sum %>% round(2), "G, ",nrow(target_table)," files."))
  urls <- target_table$url
  multi_download(urls, destfiles = file.path(database_path,urls %>% str_remove(".*/")), resume = T) %>%
    return()
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
#'   file_path = "/data/blast/downloads",
#'   database_path = "/data/blast/db",
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
  if(show_elapsed) on.exit(show_time(sec = as.numeric(Sys.time() - start_time , units = "secs") , mode =1), add = T)
  show_time("Starting to decompress the files...")
  message("[IMPORTANT]: Please ensure that the input path is a complete absolute path.")
  show_time("Unzipping data...")
  on.exit(plan(sequential), add = T)
  if (is.null(threads)) {
    plan(sequential)
  } else {
    plan(multisession, workers = threads)
  }

  file_l <- list.files(file_path, full.names = TRUE, pattern = "tar.gz")

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
#' hits, rather than a phylogenetic tree–based LCA.
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
#' @importFrom dplyr sym syms summarise filter group_by select mutate
#'
#' @examples
#' \dontrun{
#' lca_table <- find_lca(
#'   df = blast_results,
#'   seqs_ID = "qseqid",
#'   tax_levels = c("kingdom","phylum","class","order","family","genus","species")
#' )
#' }
#' @section Citations: Wang, Y., Korneliussen, T. S., Holman, L. E., Manica, A., & Pedersen, M. W. (2022). ngsLCA—A toolkit for fast and flexible lowest common ancestor inference and taxonomic profiling of metagenomic data. Methods in Ecology and Evolution, 13, 2699–2708. \url{https://doi.org/10.1111/2041-210X.14006}
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
               rank != tax_levels[length(tax_levels)] ~ paste(scientific_name, "sp.", 1:length(cur_group_rows())),.default = scientific_name
             ) )

  summ_df %>% return()

}
