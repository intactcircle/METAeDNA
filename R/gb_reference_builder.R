# GenBank reference building ---------------------------------------------------------------

if (getRversion() >= "2.15.1") utils::globalVariables(".data")

# Download NCBI BLAST databases ------------------------------------------------

#' Download the NCBI nt BLAST database
#'
#' Download all preformatted version 5 `nt` database volumes, validate them
#' once against NCBI MD5 files, download failed volumes again once, and remove
#' the downloaded MD5 files after validation.
#'
#' @param database_path Directory used for downloaded archives and metadata.
#' @param resume Resume incomplete archive files using HTTP Range requests.
#' @param multiplex Enable HTTP/2 multiplexing when supported by libcurl and
#'   the server.
#' @param proxy Optional HTTP, HTTPS, or SOCKS proxy URL, for example
#'   `"http://127.0.0.1:10808"` or `"socks5h://127.0.0.1:10808"`. When omitted,
#'   libcurl can use the standard `https_proxy` environment variable.
#' @param proxy_userpwd Optional proxy credentials in `"username:password"`
#'   form. Prefer configuring credentials outside scripts when possible.
#' @param metadata_cache_days Reuse the local NCBI directory file table for
#'   this many days.
#' @param md5_threads Number of parallel processes used for MD5 validation.
#' @param show_elapsed Report elapsed time.
#'
#' The selected file table is saved as `database_file_info.csv`.
#'
#' @return A data frame containing the selected files and MD5 results.
#'
#' @examples
#' \dontrun{
#' download_blast_db(database_path = "blastdb")
#' }
#' @export
download_blast_db <- function(database_path = ".",
                              resume = TRUE,
                              multiplex = FALSE,
                              proxy = NULL,
                              proxy_userpwd = NULL,
                              metadata_cache_days = 7,
                              md5_threads = 8L,
                              show_elapsed = TRUE) {
  dir.create(database_path, recursive = TRUE, showWarnings = FALSE)
  database_path <- normalizePath(database_path, mustWork = TRUE)
  started <- Sys.time()
  if (show_elapsed) {
    on.exit(show_time(secs = as.numeric(difftime(Sys.time(), started, units = "secs")), mode = 1), add = TRUE)
  }
  file_info_path <- file.path(database_path, "database_file_info.csv")
  file_info_current <- file.exists(file_info_path) &&
    as.numeric(difftime(Sys.time(), file.info(file_info_path)$mtime, units = "days")) <=
      metadata_cache_days
  if (file_info_current) {
    target_table <- suppressMessages(readr::read_csv(file_info_path, show_col_types = FALSE))
    if (!"website_size_g" %in% names(target_table) && "size" %in% names(target_table)) {
      target_table$website_size_g <- target_table$size
    }
    target_table$size <- NULL
    file_info_current <- all(c("file_name", "date", "website_size_g", "url") %in%
                               names(target_table))
  }
  if (file_info_current) {
    show_time("Using the local NCBI nt file table cache...")
  } else {
    base_url <- "https://ftp.ncbi.nlm.nih.gov/blast/db/"
    show_time("Reading the NCBI BLAST database directory...")
    html <- tryCatch(readLines(base_url, warn = FALSE), error = function(e) "")
    file_info_table <- data.frame(
      file_name = html %>% stringr::str_remove("[<]a.href=\"") %>%
        stringr::str_remove("\".*"),
      meta_data = html %>% stringr::str_remove(".*/a[>]") %>%
        stringr::str_trim() %>% stringr::str_squish()
    ) %>%
      dplyr::mutate(
        date = stringr::str_split(.data$meta_data, " ") %>% sapply(`[`, 1),
        size = stringr::str_split(.data$meta_data, " ") %>% sapply(`[`, 3)
      ) %>%
      dplyr::filter(stringr::str_detect(.data$file_name, "tar.gz")) %>%
      dplyr::filter(!stringr::str_detect(.data$file_name, "md5$")) %>%
      dplyr::mutate(
        size_type = dplyr::case_when(
          stringr::str_detect(.data$size, "G") ~ "G",
          stringr::str_detect(.data$size, "M") ~ "M",
          stringr::str_detect(.data$size, "K") ~ "K",
          .default = "B"
        ),
        size = stringr::str_remove(.data$size, "[A-Za-z].*") %>% as.numeric(),
        size = dplyr::case_when(
          .data$size_type == "G" ~ .data$size,
          .data$size_type == "M" ~ .data$size / 1024,
          .data$size_type == "K" ~ .data$size / (1024 * 1024),
          .default = .data$size / (1024 * 1024 * 1024)
        ),
        url = paste0(base_url, .data$file_name)
      ) %>%
      dplyr::select(-dplyr::all_of(c("size_type", "meta_data")))

    target_table <- file_info_table %>%
      dplyr::filter(stringr::str_detect(.data$file_name,
                                       "^nt[.][0-9]+[.]tar[.]gz$")) %>%
      dplyr::rename(website_size_g = "size")
    readr::write_csv(target_table, file_info_path)
  }
  show_time(paste0("Estimated file size: " %>% paste0(target_table$size %>% sum %>% round(2), "G, ",nrow(target_table)," files.")))
 #  show_time(paste0("Downloading ", nrow(target_table), " nt files."))

  archive_path <- file.path(database_path, target_table$file_name)
  checksum_path <- paste0(archive_path, ".md5")
  local_bytes <- file.info(archive_path)$size
  expected_bytes <- target_table$website_size_g * 1024^3
  show_time(paste0("Existing archive size: ",
          round(sum(local_bytes, na.rm = TRUE) / 1024^3, 2), "G."))

  pending <- is.na(local_bytes) | local_bytes < expected_bytes * 0.99
  if (any(pending)) {

    show_time(paste0("Downloading or resuming ", sum(pending), " archive(s)..."))
    download_result <- curl::multi_download(
      urls = target_table$url[pending],
      destfiles = archive_path[pending],
      resume = resume,
      progress = TRUE,
      multiplex = multiplex,
      proxy = if (is.null(proxy)) "" else proxy,
      proxyuserpwd = if (is.null(proxy_userpwd)) "" else proxy_userpwd
    )
    if (anyNA(download_result$success)) {
      stop("Download interrupted by user.", call. = FALSE)
    }
  }

  show_time(paste0("Downloading all NCBI MD5 files..."))
  checksum_result <- curl::multi_download(
    urls = paste0(target_table$url, ".md5"),
    destfiles = checksum_path,
    resume = FALSE,
    progress = FALSE,
    multiplex = multiplex,
    proxy = if (is.null(proxy)) "" else proxy,
    proxyuserpwd = if (is.null(proxy_userpwd)) "" else proxy_userpwd
  )
  if (any(!checksum_result$success | !checksum_result$status_code %in% 200:299)) {
    stop("Failed to download NCBI MD5 files.", call. = FALSE)
  }
  expected_md5 <- vapply(checksum_path, function(x) {
    tolower(sub(" .*", "", readLines(x, n = 1L, warn = FALSE)))
  }, "")

  md5_threads <- max(1L, min(as.integer(md5_threads), length(archive_path)))
  md5_cluster <- NULL
  if (md5_threads > 1L && .Platform$OS.type == "windows") {
    md5_cluster <- parallel::makeCluster(md5_threads)
    on.exit(parallel::stopCluster(md5_cluster), add = TRUE)
  }
  calculate_md5 <- function(paths) {
    if (md5_threads == 1L) return(tolower(unname(tools::md5sum(paths))))
    if (.Platform$OS.type == "windows") {
      return(tolower(unlist(parallel::parLapply(
        md5_cluster, paths, function(path) unname(tools::md5sum(path))
      ), use.names = FALSE)))
    }
    tolower(unlist(parallel::mclapply(
      paths, function(path) unname(tools::md5sum(path)),
      mc.cores = min(md5_threads, length(paths))
    ), use.names = FALSE))
  }

  md5_started <- Sys.time()
  show_time(paste0("Checking MD5 for ", length(archive_path), " files using ",
          md5_threads, " process(es)..."))
  actual_md5 <- calculate_md5(archive_path)
  show_time(paste0("MD5 check finished in ",
          round(as.numeric(difftime(Sys.time(), md5_started, units = "mins")), 2),
          " minutes."))
  incorrect <- is.na(actual_md5) | actual_md5 != expected_md5

  if (any(incorrect)) {
    show_time(paste0("Removing and re-downloading ", sum(incorrect),
            " file(s) with incorrect MD5..."))
    unlink(archive_path[incorrect])
    curl::multi_download(
      urls = target_table$url[incorrect],
      destfiles = archive_path[incorrect],
      resume = FALSE,
      progress = TRUE,
      multiplex = multiplex,
      proxy = if (is.null(proxy)) "" else proxy,
      proxyuserpwd = if (is.null(proxy_userpwd)) "" else proxy_userpwd
    )
    show_time(paste0("Rechecking MD5 for ", sum(incorrect), " downloaded file(s)..."))
    actual_md5[incorrect] <- calculate_md5(archive_path[incorrect])
  }

  target_table$website_md5 <- expected_md5
  target_table$actual_md5 <- actual_md5
  target_table$md5_ok <- !is.na(actual_md5) & actual_md5 == expected_md5
  target_table$actual_size_g <- file.info(archive_path)$size / 1024^3
  verified_size <- sum(file.info(archive_path[target_table$md5_ok])$size,
                       na.rm = TRUE) / 1024^3
  show_time(paste0("MD5 verified: ", sum(target_table$md5_ok), "/",
          nrow(target_table), " files, ", round(verified_size, 2), "G on disk."))
  readr::write_csv(target_table, file_info_path)
  md5_file_count <- sum(file.exists(checksum_path))
  unlink(checksum_path)
  show_time(paste0("Removed ", md5_file_count - sum(file.exists(checksum_path)),
                   " MD5 file(s)."))
  if (any(!target_table$md5_ok)) {
    show_time(paste0(sum(!target_table$md5_ok), " file(s) still failed MD5 validation.", call. = FALSE))
  }
  target_table
}

#' Extract accession numbers from sequence identifiers
#'
#' Remove common NCBI GI prefixes, pipe-delimited annotations, version-adjacent
#' descriptions, and text following the first space. The resulting accession
#' values are used to detect sequences representing the same database record.
#'
#' @param x Character vector containing FASTA or BLAST sequence identifiers.
#'
#' @return A character vector of cleaned accession numbers.
#'
#' @noRd
parse_gb_accession <- function(x) {
  x <- stringr::str_remove(x, "^gi[|][^|]+[|][a-z]+[|]")
  x <- stringr::str_remove(x, "[|].*")
  x <- stringr::str_remove(x, " .*")
  x
}

#' Clean and deduplicate reference sequences
#'
#' Extract accession numbers from sequence names, retain the longest sequence
#' for each duplicated accession, normalize whitespace in FASTA headers, and
#' optionally remove sequences outside the requested length range.
#'
#' @param seqs A `Biostrings::DNAStringSet` containing reference sequences.
#' @param min_length Optional minimum sequence length.
#' @param max_length Optional maximum sequence length.
#'
#' @return A cleaned `Biostrings::DNAStringSet`.
#'
#' @noRd
clean_gb_reference_sequences <- function(seqs, min_length = NULL, max_length = NULL) {
  if (length(seqs) == 0) return(seqs)
  acc <- parse_gb_accession(names(seqs))
  ord <- order(acc, Biostrings::width(seqs), decreasing = TRUE)
  seqs <- seqs[ord]
  acc <- acc[ord]
  seqs <- seqs[!duplicated(acc)]
  names(seqs) <- names(seqs) %>%
    stringr::str_replace_all("[\r\n\t]", " ") %>%
    stringr::str_squish()
  w <- Biostrings::width(seqs)
  keep <- rep(TRUE, length(seqs))
  if (!is.null(min_length)) keep <- keep & w >= min_length
  if (!is.null(max_length)) keep <- keep & w <= max_length
  seqs[keep]
}

#' Retrieve nucleotide sequences from NCBI Entrez
#'
#' Search NCBI nuccore with an organism, gene/marker name, and optional extra
#' query terms, then fetch matching FASTA records in batches. This is the
#' online half of \code{build_gb_reference()} and can also be used directly when a local
#' BLAST database is not available.
#'
#' @param taxon Target taxon used as \code{[Organism]}.
#' @param gene Gene or marker name searched in gene, title, and all fields.
#' @param extra_terms Additional raw Entrez search terms.
#' @param date_from,date_to Optional publication-date range in \code{YYYY/MM/DD}.
#' @param max_records Maximum number of records to retrieve.
#' @param batch_size Number of records fetched per API request.
#' @param email NCBI contact email. Strongly recommended by NCBI.
#' @param api_key Optional NCBI API key.
#' @param threads Retained for compatibility. NCBI requests are sent
#'   sequentially to respect E-utilities rate limits.
#' @param min_length,max_length Optional sequence length filters.
#' @param output_fasta Optional FASTA output path.
#'
#' @return A \code{DNAStringSet} of fetched and de-duplicated sequences.
#' @importFrom Biostrings DNAStringSet readDNAStringSet writeXStringSet width
#' @importFrom stringr str_remove str_replace_all str_squish
#' @export
fetch_gb_reference_online <- function(taxon = NULL, gene = NULL, extra_terms = NULL,
                              date_from = NULL, date_to = NULL,
                              max_records = 500, batch_size = 100,
                              email = NULL, api_key = NULL, threads = 1,
                              min_length = NULL, max_length = NULL,
                              output_fasta = NULL) {
  terms <- character()
  if (!is.null(taxon)) terms <- c(terms, sprintf("%s[Organism]", taxon))
  if (!is.null(gene)) terms <- c(terms, sprintf("(%s[Gene] OR %s[Title] OR %s[All Fields])", gene, gene, gene))
  if (!is.null(date_from) || !is.null(date_to)) {
    if (is.null(date_from)) date_from <- "1900/01/01"
    if (is.null(date_to)) date_to <- format(Sys.Date(), "%Y/%m/%d")
    terms <- c(terms, sprintf("%s:%s[PDAT]", date_from, date_to))
  }
  term <- paste(c(terms, extra_terms), collapse = " AND ")
  show_time(paste0("NCBI Entrez search: ", term))

  search_params <- c(
    db = "nuccore", term = term, retmode = "json",
    retmax = if (is.infinite(max_records)) "0" else as.character(max_records),
    tool = "METAeDNA"
  )
  if (!is.null(email)) search_params <- c(search_params, email = email)
  if (!is.null(api_key)) search_params <- c(search_params, api_key = api_key)
  search_url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
    paste(names(search_params), vapply(search_params, utils::URLencode, "", reserved = TRUE),
          sep = "=", collapse = "&")
  )
  search_response <- curl::curl_fetch_memory(search_url)
  search <- jsonlite::fromJSON(rawToChar(search_response$content))$esearchresult
  if (is.infinite(max_records)) {
    search_params["retmax"] <- search$count
    search_url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
      paste(names(search_params), vapply(search_params, utils::URLencode, "", reserved = TRUE),
            sep = "=", collapse = "&")
    )
    search <- jsonlite::fromJSON(
      rawToChar(curl::curl_fetch_memory(search_url)$content)
    )$esearchresult
  }
  ids <- search$idlist
  if (length(ids) == 0) {
    warning("No NCBI records matched the query.", call. = FALSE)
    return(Biostrings::DNAStringSet())
  }

  chunks <- split(ids, ceiling(seq_along(ids) / batch_size))
  show_time(paste0("Fetching ", length(ids), " NCBI records in ", length(chunks), " batch(es)."))
  fasta <- character(length(chunks))
  for (i in seq_along(chunks)) {
    fetch_params <- c(
      db = "nuccore", id = paste(chunks[[i]], collapse = ","),
      rettype = "fasta", retmode = "text", tool = "METAeDNA"
    )
    if (!is.null(email)) fetch_params <- c(fetch_params, email = email)
    if (!is.null(api_key)) fetch_params <- c(fetch_params, api_key = api_key)
    fetch_url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
      paste(names(fetch_params), vapply(fetch_params, utils::URLencode, "", reserved = TRUE),
            sep = "=", collapse = "&")
    )
    fasta[i] <- rawToChar(curl::curl_fetch_memory(fetch_url)$content)
    if (i < length(chunks)) Sys.sleep(if (is.null(api_key)) 0.34 else 0.11)
  }

  path <- tempfile(fileext = ".fasta")
  on.exit(unlink(path), add = TRUE)
  writeLines(fasta, path, useBytes = TRUE)
  seqs <- Biostrings::readDNAStringSet(path)
  seqs <- clean_gb_reference_sequences(seqs, min_length, max_length)

  if (!is.null(output_fasta)) {
    dir.create(dirname(output_fasta), recursive = TRUE, showWarnings = FALSE)
    Biostrings::writeXStringSet(seqs, output_fasta)
  }
  seqs
}

#' Build reference sequences from a local BLAST database
#'
#' Use seed sequences to search a local BLAST+ nucleotide database, retain
#' high-quality novel hits, convert the best candidates into the next query set,
#' and repeat until no new subject accessions are found or \code{max_rounds} is
#' reached.
#'
#' @param seed_fasta FASTA file containing initial query sequences.
#' @param database Local BLAST database prefix, or the data frame returned by
#'   `download_blast_db()` when one database was requested.
#' @param output_dir Directory for round-wise BLAST tables and final FASTA.
#' @param blast_exe_path Optional BLAST+ executable directory.
#' @param max_rounds Maximum iterative BLAST rounds.
#' @param max_target_seqs Maximum target sequences per query per round.
#' @param round_seed_limit Maximum number of novel candidate sequences carried
#'   forward into the next BLAST round. Candidates are ranked by an adaptive
#'   score combining bitscore, identity, and query coverage.
#' @param perc_identity Minimum BLAST percent identity.
#' @param qcov_hsp_perc Minimum BLAST query coverage.
#' @param evalue BLAST e-value threshold.
#' @param threads BLAST threads.
#' @param min_length,max_length Optional sequence length filters.
#' @param sql_path Optional local NCBI taxonomy SQLite database.
#' @param taxon Optional taxon whitelist applied through \code{sql_path}.
#' @param taxidlist Optional file containing one NCBI taxid per line. It is
#'   passed to BLAST with `-taxidlist` to restrict the local search.
#'
#' @return A list with \code{sequences}, \code{hits}, \code{rounds}, and
#'   \code{output_fasta}.
#'
#' @examples
#' \dontrun{
#' mined <- mine_gb_reference_local(
#'   seed_fasta = "reference/teleo_12S_seeds.fasta",
#'   database = "blast_database/nt",
#'   output_dir = "reference/teleo_12S_local",
#'   max_rounds = 3,
#'   max_target_seqs = 100,
#'   round_seed_limit = 200,
#'   perc_identity = 85,
#'   qcov_hsp_perc = 70,
#'   evalue = "1e-20",
#'   threads = 8,
#'   min_length = 40,
#'   max_length = 250
#' )
#'
#' mined$output_fasta
#' mined$rounds
#' head(mined$hits)
#' }
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet width
#' @importFrom dplyr arrange desc distinct filter bind_rows
#' @importFrom readr read_tsv write_tsv
#' @export
mine_gb_reference_local <- function(seed_fasta, database, output_dir = "build_gb_reference",
                             blast_exe_path = NULL, max_rounds = 3,
                             max_target_seqs = 500, round_seed_limit = 200,
                             perc_identity = 85,
                             qcov_hsp_perc = 70, evalue = "1e-20",
                             threads = 1, min_length = NULL, max_length = NULL,
                             sql_path = NULL, taxon = NULL, taxidlist = NULL) {
  if (is.data.frame(database)) database <- unique(database$database_prefix)[1]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  current <- Biostrings::readDNAStringSet(seed_fasta)
  all_sequences <- current
  all_hits <- list()
  seen <- parse_gb_accession(names(current))

  for (round_id in seq_len(max_rounds)) {
    show_time(paste0("Local BLAST reference building round ", round_id, "/", max_rounds))
    query_file <- file.path(output_dir, sprintf("round_%02d_query.fasta", round_id))
    out_file <- file.path(output_dir, sprintf("round_%02d_blast.tsv", round_id))
    Biostrings::writeXStringSet(current, query_file)

    blast_cmd(
      base_cmd = "blastn",
      args = c(
        "-query", query_file, "-db", database,
        "-outfmt", shQuote("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxids sseq"),
        "-max_target_seqs", max_target_seqs, "-out", out_file,
        "-perc_identity", perc_identity, "-qcov_hsp_perc", qcov_hsp_perc,
        "-evalue", evalue, "-num_threads", threads,
        if (!is.null(taxidlist)) c("-taxidlist", taxidlist)
      ),
      blast_exe_path = blast_exe_path
    )

    if (!file.exists(out_file) || file.info(out_file)$size == 0) break
    columns <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                 "qcovs", "staxids", "sseq")
    hits <- readr::read_tsv(out_file, col_names = columns, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::filter(.data$pident >= perc_identity, .data$qcovs >= qcov_hsp_perc)
    all_hits[[round_id]] <- hits
    hits <- hits %>%
      dplyr::mutate(
        accession = parse_gb_accession(.data$sseqid),
        score = .data$bitscore * .data$pident * .data$qcovs
      ) %>%
      dplyr::filter(!.data$accession %in% seen) %>%
      dplyr::arrange(dplyr::desc(.data$score)) %>%
      dplyr::distinct(.data$accession, .keep_all = TRUE)
    if (!is.null(round_seed_limit)) hits <- utils::head(hits, round_seed_limit)
    if (!nrow(hits)) {
      show_time("No new accessions found; stopping iteration.")
      break
    }

    round_sequences <- Biostrings::DNAStringSet(hits$sseq)
    names(round_sequences) <- paste(hits$sseqid, "taxid", hits$staxids, sep = "|")
    round_sequences <- clean_gb_reference_sequences(round_sequences, min_length, max_length)

    if (!is.null(sql_path) && !is.null(taxon) && length(round_sequences)) {
      ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
      taxid <- get_taxid(sql_path, parse_gb_accession(names(round_sequences)), show_elapsed = FALSE)
      taxonomy <- get_taxonomy(sql_path, unique(stats::na.omit(taxid$taxid)),
                               rank_name = ranks, show_elapsed = FALSE)
      taxonomy <- dplyr::left_join(taxid, taxonomy, by = dplyr::join_by(taxid == query_taxid))
      keep <- apply(taxonomy[, intersect(ranks, names(taxonomy)), drop = FALSE], 1, function(x) taxon %in% x)
      round_sequences <- round_sequences[keep]
    }

    acc <- parse_gb_accession(names(round_sequences))
    new_idx <- !(acc %in% seen)
    if (!any(new_idx)) {
      show_time("No new accessions found; stopping iteration.")
      break
    }
    current <- round_sequences[new_idx]
    all_sequences <- c(all_sequences, current)
    seen <- unique(c(seen, acc[new_idx]))
  }

  all_sequences <- clean_gb_reference_sequences(all_sequences, min_length, max_length)
  output_fasta <- file.path(output_dir, "build_gb_reference_sequences.fasta")
  Biostrings::writeXStringSet(all_sequences, output_fasta)

  hits_tbl <- if (length(all_hits) > 0) dplyr::bind_rows(all_hits, .id = "round") else data.frame()
  if (nrow(hits_tbl) > 0) readr::write_tsv(hits_tbl, file.path(output_dir, "build_gb_reference_hits.tsv"))

  list(
    sequences = all_sequences,
    hits = hits_tbl,
    rounds = length(all_hits),
    output_fasta = output_fasta
  )
}

#' Build GenBank reference sequences online or locally
#'
#' Online mode retrieves records from NCBI E-utilities using a taxon and marker
#' query. Local mode starts from a seed FASTA and expands candidate sequences by
#' iterative BLAST against a local nucleotide database. Both modes clean,
#' de-duplicate, filter, and write a reusable FASTA collection.
#'
#' @param mode One of \code{"online"} or \code{"local"}.
#' @param taxon,gene,extra_terms,date_from,date_to Entrez search settings for
#'   online mode. The taxon can also be used for validation in local mode.
#' @param seed_fasta Initial FASTA for local mode.
#' @param database Local BLAST database prefix for local mode, or
#'   the result returned by `download_blast_db()` for one database.
#' @param output_dir Output directory.
#' @param email,api_key NCBI contact/API settings.
#' @param max_records,batch_size Online retrieval controls.
#' @param max_rounds,max_target_seqs,round_seed_limit,perc_identity,qcov_hsp_perc,evalue
#'   Local iterative BLAST controls.
#' @param threads Number of API fetch workers or BLAST threads.
#' @param min_length,max_length Optional sequence length filters.
#' @param sql_path Optional local NCBI taxonomy SQLite database for taxon
#'   validation in local mode.
#' @param blast_exe_path Optional BLAST+ executable directory.
#'
#' @return A list containing final \code{sequences}, output file paths, and
#'   mode-specific metadata.
#' @examples
#' \dontrun{
#' res <- build_gb_reference(
#'   mode = "online",
#'   taxon = "Saxifragales",
#'   gene = "ITS",
#'   email = "your.name@example.org",
#'   max_records = 200,
#'   threads = 2
#' )
#'
#' res2 <- build_gb_reference(
#'   mode = "local",
#'   seed_fasta = "example_data/seed_sequences.fasta",
#'   database = "reference_databases/blast/nt",
#'   output_dir = "analysis_results/ITS_Saxifragales",
#'   threads = 8
#' )
#'
#' # The download result can be passed directly after decompression.
#' db <- download_blast_db("nt_euk", "reference_databases/blast",
#'                         decompress = TRUE, keep_archives = FALSE)
#' res3 <- build_gb_reference(
#'   mode = "local", seed_fasta = "seed.fasta", database = db,
#'   output_dir = "analysis_results/reference", threads = 8
#' )
#' }
#' @export
build_gb_reference <- function(mode = c("online", "local"),
                       taxon = NULL, gene = NULL, extra_terms = NULL,
                       date_from = NULL, date_to = NULL,
                       seed_fasta = NULL, database = NULL,
                       output_dir = "build_gb_reference",
                       email = NULL, api_key = NULL,
                       max_records = 500, batch_size = 100,
                       max_rounds = 3, max_target_seqs = 50,
                       round_seed_limit = 200,
                       perc_identity = 85, qcov_hsp_perc = 70,
                       evalue = "1e-20", threads = 1,
                       min_length = NULL, max_length = NULL,
                       sql_path = NULL, blast_exe_path = NULL) {
  mode <- match.arg(mode)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  start_time <- Sys.time()
  on.exit(show_time(secs = as.numeric(Sys.time() - start_time, units = "secs"), mode = 1), add = TRUE)

  if (mode == "online") {
    output_fasta <- file.path(output_dir, "build_gb_reference_online.fasta")
    seqs <- fetch_gb_reference_online(
      taxon = taxon, gene = gene, extra_terms = extra_terms,
      date_from = date_from, date_to = date_to,
      max_records = max_records, batch_size = batch_size,
      email = email, api_key = api_key, threads = threads,
      min_length = min_length, max_length = max_length,
      output_fasta = output_fasta
    )
    return(list(mode = mode, sequences = seqs, output_fasta = output_fasta))
  }

  local <- mine_gb_reference_local(
    seed_fasta = seed_fasta,
    database = database,
    output_dir = output_dir,
    blast_exe_path = blast_exe_path,
    max_rounds = max_rounds,
    max_target_seqs = max_target_seqs,
    round_seed_limit = round_seed_limit,
    perc_identity = perc_identity,
    qcov_hsp_perc = qcov_hsp_perc,
    evalue = evalue,
    threads = threads,
    min_length = min_length,
    max_length = max_length,
    sql_path = sql_path,
    taxon = taxon
  )
  c(list(mode = mode), local)
}

#' Build a local nucleotide reference database
#'
#' Build a marker-specific nucleotide database from one of three sources: the
#' result of `col_ncbi_crosswalk()`, an NCBI taxonomic node supplied as a name or
#' taxid, or reference sequences expanded iteratively against a local NCBI BLAST
#' database. Online modes download matching NCBI Nucleotide records. Seed mode
#' behaves like NCBIminer and can restrict BLAST to all descendants of selected
#' NCBI nodes using the local taxonomy SQLite database.
#'
#' Common marker annotations are translated automatically. Examples include
#' `trnL` (chloroplast trnL and trnL-trnF), `rbcL`, `matK`, `ITS` (ITS1/ITS2),
#' `COI`/`CO1`, `12S`, `16S`/small-subunit rRNA,
#' `18S`/eukaryotic small-subunit rRNA, and `28S`/large-subunit rRNA. Supply
#' `marker_query` to replace the built-in expression for another annotation.
#'
#' @param crosswalk Optional table returned by `col_ncbi_crosswalk()`.
#' @param marker Optional marker annotation, for example `"trnL"`, `"ITS"`, or
#'   `"rbcL"`. When omitted, all nucleotide records below the selected taxon are
#'   eligible.
#' @param taxon Optional NCBI taxon name or vector of names. Descendants are
#'   included.
#' @param taxid Optional NCBI taxid or vector of taxids. Descendants are
#'   included.
#' @param seed_fasta Optional FASTA containing one or more reference sequences.
#'   When supplied, iterative local BLAST mode is used.
#' @param search_database Local NCBI nucleotide BLAST database prefix used in
#'   seed mode, for example `"blastdb/nt_euk"`.
#' @param ncbi_db Local NCBI taxonomy SQLite database. It is required when seed
#'   mode is restricted by `taxon`, `taxid`, or `crosswalk`.
#' @param marker_query Optional raw NCBI Entrez expression overriding the
#'   built-in marker expression.
#' @param output_dir Output directory.
#' @param max_records Maximum online records. `Inf` downloads every record
#'   reported by the query.
#' @param batch_size Number of NCBI records fetched per request.
#' @param max_rounds,max_target_seqs,round_seed_limit Iterative BLAST controls.
#' @param perc_identity,qcov_hsp_perc,evalue BLAST match thresholds.
#' @param min_length,max_length Optional final sequence-length limits.
#' @param threads Number of local BLAST threads.
#' @param email,api_key NCBI E-utilities identification and API settings.
#' @param blast_exe_path Optional BLAST+ executable directory.
#' @param makeblastdb If `TRUE`, format the final FASTA as a nucleotide BLAST
#'   database.
#' @param database_name Basename used for the final FASTA and BLAST database.
#'
#' @return A list containing `mode`, `query`, `sequences`, `fasta`,
#'   `blast_database`, and source-specific results.
#' @export
#'
#' @examples
#' \dontrun{
#' # Download trnL records for taxa resolved by COL and NCBI.
#' build_nucleotide_database(crosswalk = matched_names, marker = "trnL")
#'
#' # Download all rbcL records below a taxonomic node.
#' build_nucleotide_database(
#'   taxon = "Saxifragales", marker = "rbcL", max_records = Inf
#' )
#'
#' # NCBIminer-like expansion restricted to a node.
#' build_nucleotide_database(
#'   seed_fasta = "seed.fasta", search_database = "blastdb/nt_euk",
#'   taxid = 41944, ncbi_db = "ncbi_taxonomy.sqlite", threads = 16
#' )
#' }
build_nucleotide_database <- function(crosswalk = NULL,
                                      marker = NULL,
                                      taxon = NULL,
                                      taxid = NULL,
                                      seed_fasta = NULL,
                                      search_database = NULL,
                                      ncbi_db = NULL,
                                      marker_query = NULL,
                                      output_dir = "nucleotide_database",
                                      max_records = Inf,
                                      batch_size = 200,
                                      max_rounds = 3,
                                      max_target_seqs = 100,
                                      round_seed_limit = 500,
                                      perc_identity = 85,
                                      qcov_hsp_perc = 70,
                                      evalue = "1e-20",
                                      min_length = NULL,
                                      max_length = NULL,
                                      threads = 1,
                                      email = NULL,
                                      api_key = NULL,
                                      blast_exe_path = NULL,
                                      makeblastdb = TRUE,
                                      database_name = "reference") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  crosswalk_taxids <- if (!is.null(crosswalk) && "ncbi_taxid" %in% names(crosswalk)) {
    unique(stats::na.omit(as.integer(crosswalk$ncbi_taxid)))
  } else integer()
  selected_taxids <- unique(c(as.integer(taxid), crosswalk_taxids))
  selected_taxids <- selected_taxids[!is.na(selected_taxids)]

  if (!is.null(seed_fasta)) {
    taxidlist <- NULL
    restriction <- length(selected_taxids) || length(taxon)
    if (restriction) {
      con <- DBI::dbConnect(RSQLite::SQLite(), ncbi_db)
      on.exit(DBI::dbDisconnect(con), add = TRUE)
      if (length(taxon)) {
        name_table <- if ("name_aliases" %in% DBI::dbListTables(con)) "name_aliases" else "names"
        name_ids <- DBI::dbGetQuery(
          con,
          sprintf(
            "SELECT DISTINCT tax_id FROM %s WHERE name IN (%s)",
            name_table,
            paste(DBI::dbQuoteString(con, taxon), collapse = ",")
          )
        )$tax_id
        selected_taxids <- unique(c(selected_taxids, name_ids))
      }
      roots <- paste(DBI::dbQuoteLiteral(con, selected_taxids), collapse = ",")
      descendants <- DBI::dbGetQuery(
        con,
        paste0(
          "WITH RECURSIVE descendants(tax_id) AS (",
          " SELECT tax_id FROM nodes WHERE tax_id IN (", roots, ")",
          " UNION SELECT n.tax_id FROM nodes n JOIN descendants d",
          " ON n.parent_tax_id = d.tax_id WHERE n.tax_id <> n.parent_tax_id",
          ") SELECT DISTINCT tax_id FROM descendants"
        )
      )$tax_id
      taxidlist <- file.path(output_dir, "restricted_taxids.txt")
      writeLines(as.character(descendants), taxidlist)
    }

    source_result <- mine_gb_reference_local(
      seed_fasta = seed_fasta,
      database = search_database,
      output_dir = file.path(output_dir, "mining"),
      blast_exe_path = blast_exe_path,
      max_rounds = max_rounds,
      max_target_seqs = max_target_seqs,
      round_seed_limit = round_seed_limit,
      perc_identity = perc_identity,
      qcov_hsp_perc = qcov_hsp_perc,
      evalue = evalue,
      threads = threads,
      min_length = min_length,
      max_length = max_length,
      taxidlist = taxidlist
    )
    sequences <- source_result$sequences
    query <- if (is.null(taxidlist)) "local BLAST" else paste("local BLAST restricted by", taxidlist)
    mode <- "seed"
  } else {
    marker_key <- tolower(gsub("[^a-z0-9]", "", marker))
    marker_key[marker_key %in% c("co1", "cox1")] <- "coi"
    marker_queries <- c(
      trnl = '(trnL[Gene] OR "trnL-trnF"[Title] OR "tRNA-Leu"[Title])',
      rbcl = '(rbcL[Gene] OR "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"[Title])',
      matk = '(matK[Gene] OR "maturase K"[Title])',
      its = '(ITS1[Title] OR ITS2[Title] OR "internal transcribed spacer"[Title] OR ITS)',
      coi = '(COI[Title] OR CO1[Title] OR COX1[Gene] OR "cytochrome c oxidase subunit I"[Title])',
      `12s` = '("12S ribosomal RNA"[Title] OR 12S[Title])',
      `16s` = '("16S ribosomal RNA"[Title] OR "16S rRNA"[Title] OR 16S[Title] OR "small subunit ribosomal RNA"[Title] OR "small subunit rRNA"[Title] OR SSU[Title])',
      `18s` = '("18S ribosomal RNA"[Title] OR "18S rRNA"[Title] OR 18S[Title] OR "eukaryotic small subunit ribosomal RNA"[Title] OR "nuclear small subunit ribosomal RNA"[Title] OR "small subunit ribosomal RNA"[Title] OR "small subunit rRNA"[Title] OR "SSU rRNA"[Title] OR SSU[Title])',
      `28s` = '("28S ribosomal RNA"[Title] OR "28S rRNA"[Title] OR 28S[Title] OR "large subunit ribosomal RNA"[Title] OR "large subunit rRNA"[Title] OR LSU[Title])'
    )
    if (is.null(marker_query) && length(marker_key)) {
      marker_query <- unname(marker_queries[marker_key])
      if (is.na(marker_query)) {
        marker_query <- sprintf('(%s[Gene] OR "%s"[Title] OR "%s"[All Fields])',
                                marker, marker, marker)
      }
    }

    organism_terms <- c(
      if (length(selected_taxids)) paste0("txid", selected_taxids, "[Organism:exp]"),
      if (length(taxon)) paste0('"', taxon, '"[Organism:exp]')
    )
    if (!length(organism_terms) && !is.null(crosswalk)) {
      name_column <- if ("col_accepted_name" %in% names(crosswalk)) "col_accepted_name" else "standardized_name"
      organism_terms <- paste0('"', unique(stats::na.omit(crosswalk[[name_column]])), '"[Organism]')
    }
    query_parts <- c(
      if (length(organism_terms)) paste0("(", paste(organism_terms, collapse = " OR "), ")"),
      marker_query
    )
    query <- paste(query_parts, collapse = " AND ")
    source_result <- fetch_gb_reference_online(
      extra_terms = query,
      max_records = max_records,
      batch_size = batch_size,
      email = email,
      api_key = api_key,
      min_length = min_length,
      max_length = max_length
    )
    sequences <- source_result
    mode <- "online"
  }

  sequences <- clean_gb_reference_sequences(sequences, min_length, max_length)
  fasta <- file.path(output_dir, paste0(database_name, ".fasta"))
  Biostrings::writeXStringSet(sequences, fasta)

  blast_database <- NULL
  if (isTRUE(makeblastdb)) {
    blast_database <- file.path(output_dir, database_name)
    blast_cmd(
      base_cmd = "makeblastdb",
      args = c("-in", fasta, "-dbtype", "nucl", "-out", blast_database,
               "-title", database_name),
      blast_exe_path = blast_exe_path
    )
  }

  list(
    mode = mode,
    query = query,
    sequences = sequences,
    fasta = fasta,
    blast_database = blast_database,
    source = source_result
  )
}
