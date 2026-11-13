# GenBank reference building ---------------------------------------------------------------

if (getRversion() >= "2.15.1") utils::globalVariables(".data")

#' Common NCBI marker queries
#'
#' Entrez expressions for retrieving candidate nucleotide records. These are
#' supplied through `ncbi_query` in `fetch_gb_reference_online()` or
#' `build_gb_reference()`. They narrow the NCBI search;
#' use `marker_regex` to check record titles afterward.
#'
#' @format A named list.
#' @export
marker_entrez_query <- list(
  trnL = '(trnL[Gene] OR trnL[Title] OR "tRNA-Leu"[Title])',
  rbcL = '(rbcL[Gene] OR "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"[Title])',
  matK = '(matK[Gene] OR "maturase K"[Title])',
  ITS = '(ITS1[Title] OR ITS2[Title] OR "internal transcribed spacer"[Title] OR ITS)',
  COI = '(COI[Title] OR CO1[Title] OR COX1[Gene] OR "cytochrome c oxidase subunit I"[Title])',
  `12S` = '("12S ribosomal RNA"[Title] OR 12S[Title])',
  `16S` = '("16S ribosomal RNA"[Title] OR "16S rRNA"[Title] OR 16S[Title])',
  `18S` = '("18S ribosomal RNA"[Title] OR "18S rRNA"[Title] OR 18S[Title])',
  `28S` = '("28S ribosomal RNA"[Title] OR "28S rRNA"[Title] OR 28S[Title])'
)

#' Common marker title regular expressions
#'
#' Regular expressions for checking NCBI FASTA or BLAST record titles. Supply
#' one explicitly through `marker` or `regex_expression`.
#'
#' @format A named list.
#' @export
marker_regex <- list(
  rbcL = c("\\brbcL\\b", "ribulose-1,5-bisphosphate carboxylase"),
  matK = c("\\bmatK\\b", "maturase K"),
  trnL = paste0(
    "^(?!.*mitochond)(?!.*trnL-trnF\\s+intergenic\\s+spacer)",
    "(?=.*(?:\\btrnL\\b|tRNA-Leu))",
    "(?=.*(?:chloroplast|plastid|\\bintron\\b)).*$"
  ),
  COI = c("\\b(?:COI|CO1|COX1)\\b", "cytochrome c oxidase subunit I"),
  `12S` = c("\\b12S\\b", "12S ribosomal RNA"),
  `16S` = c("\\b16S\\b", "16S ribosomal RNA"),
  `18S` = c("\\b18S\\b", "18S ribosomal RNA"),
  `28S` = c("\\b28S\\b", "28S ribosomal RNA"),
  ITS = c("\\bITS(?:1|2)?\\b", "internal transcribed spacer")
)

# Download NCBI BLAST databases ------------------------------------------------

#' Download the NCBI nt BLAST database
#'
#' Download or resume all preformatted version 5 `nt` database volumes. After
#' every archive transfer finishes, validate all archives against the NCBI MD5
#' files, download failed archives again once, and remove every `.md5` file
#' from the download directory when validation succeeds.
#'
#' @param database_path Directory used for downloaded archives and metadata.
#' @param resume Resume `.part` archive files using HTTP Range requests. A
#'   resumed response is accepted only when the server confirms a byte-range
#'   response; otherwise the affected archive is downloaded from the beginning.
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
#'   File-level progress is displayed while checksums are calculated.
#'
#' @details The selected file table is saved as `database_file_info.csv`.
#' Existing archives and `.part` files are submitted as resume candidates;
#' their integrity is determined only by the final MD5 check after all archive
#' transfers have finished. Elapsed time is reported after successful completion.
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
                              md5_threads = 8L) {
  .check_path(database_path, "database_path")
  .check_flag(resume, "resume")
  .check_flag(multiplex, "multiplex")
  .check_number(metadata_cache_days, "metadata_cache_days", 0)
  .check_number(md5_threads, "md5_threads", 1, .Machine$integer.max, TRUE)
  dir.create(database_path, recursive = TRUE, showWarnings = FALSE)
  database_path <- normalizePath(database_path, mustWork = TRUE)
  started <- Sys.time()
  file_info_path <- file.path(database_path, "database_file_info.csv")
  file_info_current <- file.exists(file_info_path) &&
    as.numeric(difftime(Sys.time(), file.info(file_info_path)$mtime, units = "days")) <=
      metadata_cache_days
  if (file_info_current) {
    target_table <- tryCatch(
      suppressMessages(readr::read_csv(file_info_path, show_col_types = FALSE)),
      error = function(e) NULL
    )
    if (!is.null(target_table) &&
        !"website_size_g" %in% names(target_table) &&
        "size" %in% names(target_table)) {
      target_table$website_size_g <- target_table$size
    }
    if (!is.null(target_table) && "size" %in% names(target_table)) {
      target_table$size <- NULL
    }
    required_columns <- c("file_name", "date", "website_size_g", "url")
    file_info_current <- !is.null(target_table) &&
      all(required_columns %in% names(target_table)) &&
      nrow(target_table) > 0L
    if (file_info_current) {
      target_table$website_size_g <- suppressWarnings(
        as.numeric(target_table$website_size_g)
      )
      file_info_current <- all(
        grepl("^nt[.][0-9]+[.]tar[.]gz$", target_table$file_name) &
          is.finite(target_table$website_size_g) &
          target_table$website_size_g > 0 &
          !is.na(target_table$url) & nzchar(target_table$url)
      )
    }
  }
  if (file_info_current) {
    show_time("Using the local NCBI nt file table cache...")
  } else {
    base_url <- "https://ftp.ncbi.nlm.nih.gov/blast/db/"
    show_time("Reading the NCBI BLAST database directory...")
    request <- curl::new_handle(
      timeout = 120L, connecttimeout = 30L, followlocation = TRUE
    )
    if (!is.null(proxy)) curl::handle_setopt(request, proxy = proxy)
    if (!is.null(proxy_userpwd)) {
      curl::handle_setopt(request, proxyuserpwd = proxy_userpwd)
    }
    html <- tryCatch(
      rawToChar(curl::curl_fetch_memory(base_url, handle = request)$content),
      error = function(e) ""
    )
    html <- strsplit(html, "\n", fixed = TRUE)[[1L]]
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
    if (!nrow(target_table) ||
        any(!is.finite(target_table$website_size_g)) ||
        any(target_table$website_size_g <= 0)) {
      stop(
        "No valid NCBI nt archive information was found. The remote directory may be unavailable or its format may have changed.",
        call. = FALSE
      )
    }
    readr::write_csv(target_table, file_info_path)
  }
  show_time(paste0(
    "Estimated file size: ", round(sum(target_table$website_size_g), 2),
    "G, ", nrow(target_table), " files."
  ))
 #  show_time(paste0("Downloading ", nrow(target_table), " nt files."))

  archive_path <- file.path(database_path, target_table$file_name)
  checksum_path <- paste0(archive_path, ".md5")
  part_path <- paste0(archive_path, ".part")

  for (i in seq_along(archive_path)) {
    if (file.exists(archive_path[i])) {
      unlink(part_path[i])
      if (isTRUE(resume)) {
        if (!file.rename(archive_path[i], part_path[i])) {
          stop("Cannot prepare the partial download file: ", part_path[i],
               call. = FALSE)
        }
      } else {
        unlink(archive_path[i])
      }
    }
    if (!isTRUE(resume)) unlink(part_path[i])
  }

  show_time(paste0("Downloading or resuming all ", length(archive_path),
                   " archive(s)..."))
  had_partial <- isTRUE(resume) & file.exists(part_path)
  download_result <- curl::multi_download(
    urls = target_table$url,
    destfiles = part_path,
    resume = resume,
    progress = TRUE,
    multiplex = multiplex,
    proxy = if (is.null(proxy)) "" else proxy,
    proxyuserpwd = if (is.null(proxy_userpwd)) "" else proxy_userpwd
  )
  status_code <- download_result$status_code
  range_complete <- had_partial & !is.na(status_code) &
    status_code == 416L & file.exists(part_path)
  download_ok <- (!is.na(download_result$success) &
                    download_result$success &
                    ((!had_partial & status_code %in% 200:299) |
                       (had_partial & status_code == 206L))) |
    range_complete
  download_ok[is.na(download_ok)] <- FALSE

  unsafe_resume <- had_partial & !is.na(status_code) & status_code == 200L
  if (any(unsafe_resume)) {
    show_time(paste0("The server restarted ", sum(unsafe_resume),
                     " resumed download(s); downloading them from the beginning..."))
    unlink(part_path[unsafe_resume])
    full_result <- curl::multi_download(
      urls = target_table$url[unsafe_resume],
      destfiles = part_path[unsafe_resume],
      resume = FALSE,
      progress = TRUE,
      multiplex = multiplex,
      proxy = if (is.null(proxy)) "" else proxy,
      proxyuserpwd = if (is.null(proxy_userpwd)) "" else proxy_userpwd
    )
    full_ok <- !is.na(full_result$success) & full_result$success &
      full_result$status_code %in% 200:299
    full_ok[is.na(full_ok)] <- FALSE
    download_ok[unsafe_resume] <- full_ok
  }
  if (any(!download_ok)) {
    stop(sum(!download_ok),
         " NCBI nt archive download(s) failed. Partial files were retained for the next run.",
         call. = FALSE)
  }
  for (i in seq_along(archive_path)) {
    if (!file.rename(part_path[i], archive_path[i])) {
      stop("Cannot finalize the downloaded archive: ", archive_path[i],
           call. = FALSE)
    }
  }

  show_time("Downloading all NCBI MD5 files...")
  checksum_result <- curl::multi_download(
    urls = paste0(target_table$url, ".md5"),
    destfiles = checksum_path,
    resume = FALSE,
    progress = FALSE,
    multiplex = multiplex,
    proxy = if (is.null(proxy)) "" else proxy,
    proxyuserpwd = if (is.null(proxy_userpwd)) "" else proxy_userpwd
  )
  checksum_failed <- is.na(checksum_result$success) |
    !checksum_result$success |
    !checksum_result$status_code %in% 200:299
  if (any(checksum_failed)) {
    show_time(paste0("Retrying ", sum(checksum_failed), " failed MD5 download(s)..."))
    checksum_retry <- curl::multi_download(
      urls = paste0(target_table$url[checksum_failed], ".md5"),
      destfiles = checksum_path[checksum_failed],
      resume = FALSE,
      progress = FALSE,
      multiplex = multiplex,
      proxy = if (is.null(proxy)) "" else proxy,
      proxyuserpwd = if (is.null(proxy_userpwd)) "" else proxy_userpwd
    )
    retry_failed <- is.na(checksum_retry$success) |
      !checksum_retry$success |
      !checksum_retry$status_code %in% 200:299
    if (any(retry_failed)) {
      stop(sum(retry_failed), " NCBI MD5 file download(s) failed after retry.",
           call. = FALSE)
    }
  }
  expected_md5 <- vapply(checksum_path, function(x) {
    tolower(sub(" .*", "", readLines(x, n = 1L, warn = FALSE)))
  }, "")
  if (any(!grepl("^[0-9a-f]{32}$", expected_md5))) {
    stop("One or more downloaded NCBI MD5 files are invalid.", call. = FALSE)
  }

  md5_threads <- max(1L, min(as.integer(md5_threads), length(archive_path)))
  md5_cluster <- NULL
  if (md5_threads > 1L && .Platform$OS.type == "windows") {
    md5_cluster <- parallel::makeCluster(md5_threads)
    on.exit(parallel::stopCluster(md5_cluster), add = TRUE)
  }
  calculate_md5 <- function(paths) {
    if (!length(paths)) return(character())
    workers <- min(md5_threads, length(paths))
    batches <- split(seq_along(paths),
                     ceiling(seq_along(paths) / workers))
    result <- character(length(paths))
    progress <- utils::txtProgressBar(min = 0, max = length(paths), width = 40, style = 3)
    on.exit(close(progress), add = TRUE)

    completed <- 0L
    for (index in batches) {
      if (workers == 1L) {
        value <- tools::md5sum(paths[index])
      } else if (.Platform$OS.type == "windows") {
        value <- unlist(parallel::parLapply(
          md5_cluster, paths[index], function(path) unname(tools::md5sum(path))
        ), use.names = FALSE)
      } else {
        value <- unlist(parallel::mclapply(
          paths[index], function(path) unname(tools::md5sum(path)),
          mc.cores = min(workers, length(index))
        ), use.names = FALSE)
      }
      result[index] <- tolower(unname(value))
      completed <- completed + length(index)
      utils::setTxtProgressBar(progress, completed)
    }
    result
  }

  show_time(paste0("Running the final MD5 check for ", length(archive_path),
                   " archive(s)..."))
  actual_md5 <- calculate_md5(archive_path)
  incorrect <- is.na(actual_md5) | actual_md5 != expected_md5
  if (any(incorrect)) {
    show_time(paste0("Downloading ", sum(incorrect),
                     " file(s) again from the beginning after MD5 failure..."))
    unlink(c(archive_path[incorrect], part_path[incorrect]))
    retry_result <- curl::multi_download(
      urls = target_table$url[incorrect],
      destfiles = part_path[incorrect],
      resume = FALSE,
      progress = TRUE,
      multiplex = multiplex,
      proxy = if (is.null(proxy)) "" else proxy,
      proxyuserpwd = if (is.null(proxy_userpwd)) "" else proxy_userpwd
    )
    retry_ok <- !is.na(retry_result$success) & retry_result$success &
      retry_result$status_code %in% 200:299
    retry_ok[is.na(retry_ok)] <- FALSE
    if (any(!retry_ok)) {
      stop(sum(!retry_ok), " NCBI nt archive re-download(s) failed.",
           call. = FALSE)
    }
    for (i in which(incorrect)) {
      if (!file.rename(part_path[i], archive_path[i])) {
        stop("Cannot finalize the downloaded archive: ", archive_path[i],
             call. = FALSE)
      }
    }
    show_time(paste0("Rechecking MD5 for ", sum(incorrect),
                     " re-downloaded archive(s)..."))
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
  if (any(!target_table$md5_ok)) {
    stop(sum(!target_table$md5_ok),
         " file(s) still failed MD5 validation. MD5 files were retained.",
         call. = FALSE)
  }
  show_time("Final MD5 check passed for all archives.")
  md5_files <- list.files(database_path, pattern = "[.]md5$",
                          full.names = TRUE, ignore.case = TRUE)
  unlink(md5_files)
  remaining_md5 <- file.exists(md5_files)
  show_time(paste0("Removed ", sum(!remaining_md5), " MD5 file(s)."))
  if (any(remaining_md5)) {
    warning(sum(remaining_md5), " MD5 file(s) could not be removed.",
            call. = FALSE)
  }
  show_time(secs = as.numeric(difftime(Sys.time(), started, units = "secs")), mode = 1)
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
#' Search NCBI nuccore by taxon name or TaxID using Entrez syntax. `ncbi_query`
#' adds an NCBI search expression (for example, a marker). The optional
#' `regex_expression` then checks record titles from ESummary. After all titles
#' have been checked,
#' matching accessions are downloaded in batches as FASTA by EFetch. A record
#' excluded by the Entrez query cannot be recovered
#' by the regular expression, so choose a broad candidate query.
#' ESearch pages through candidate UIDs, which are sent directly to ESummary.
#'
#' @param taxon_names NCBI taxon names. Descendants are included.
#' @param taxids NCBI TaxIDs or a table with an `ncbi_taxid` column, such as
#'   the result of `col_ncbi_crosswalk()`. Descendants are included.
#' @param gene Gene or marker name searched in gene, title, and all fields.
#' @param ncbi_query Optional NCBI Entrez query for candidate records, such as
#'   `marker_entrez_query$trnL`. Omitting it searches all records for the
#'   selected taxa; it is not an R regular expression.
#' @param regex_expression Optional R regular expression(s) matched against
#'   NCBI record titles before FASTA download, such as `marker_regex$trnL`.
#'   Multiple expressions are combined with OR.
#' @param date_from,date_to Optional publication-date range in \code{YYYY/MM/DD}.
#' @param max_records Maximum number of Entrez candidate records to inspect.
#'   Use `Inf` to check all candidates returned by the NCBI query.
#' @param batch_size Maximum accessions per EFetch FASTA request. Titles are
#'   retrieved separately in batches of up to 10,000 XML summaries.
#' @param email NCBI contact email. Strongly recommended by NCBI.
#' @param api_key Optional NCBI API key.
#' @param n_threads Retained for compatibility. NCBI requests are sent
#'   sequentially. E-utilities requests, including retries, are limited to
#'   at most 3 per second without `api_key` or 10 per second with one.
#' @param min_length,max_length Optional sequence length filters.
#' @param output_fasta Optional FASTA output path.
#' @param taxon_count Optional number of species represented by a pre-built
#'   `ncbi_query`. It is used only in progress messages.
#'
#' @return A \code{DNAStringSet} of fetched and de-duplicated sequences. Its
#'   `taxid_map` attribute contains accession-to-TaxID pairs for `makeblastdb`.
#' @examples
#' \dontrun{
#' fetch_gb_reference_online(
#'   taxids = c(9606, 9598), ncbi_query = marker_entrez_query$COI,
#'   regex_expression = marker_regex$COI,
#'   output_fasta = "reference/COI.fasta"
#' )
#' fetch_gb_reference_online(
#'   taxon_names = c("Quercus robur", "Quercus ilex"),
#'   ncbi_query = marker_entrez_query$rbcL,
#'   regex_expression = marker_regex$rbcL
#' )
#' }
#' @importFrom Biostrings DNAStringSet readDNAStringSet writeXStringSet width
#' @importFrom stringr str_remove str_replace_all str_squish
#' @export
fetch_gb_reference_online <- function(taxon_names = NULL, taxids = NULL,
                                      gene = NULL, ncbi_query = NULL,
                                      regex_expression = NULL,
                                      date_from = NULL, date_to = NULL,
                                      max_records = 500, batch_size = 100,
                                      email = NULL, api_key = NULL, n_threads = 1,
                                      min_length = NULL, max_length = NULL,
                                      output_fasta = NULL, taxon_count = NULL) {
  if (is.data.frame(taxids)) {
    if (!"ncbi_taxid" %in% names(taxids)) {
      stop("taxids table must contain an ncbi_taxid column.", call. = FALSE)
    }
    taxids <- taxids$ncbi_taxid
  }
  taxids <- unique(trimws(as.character(taxids)))
  taxids <- taxids[!is.na(taxids) & nzchar(taxids)]
  if (length(taxids) && any(!grepl("^[1-9][0-9]*$", taxids))) {
    stop("taxids must contain positive NCBI TaxIDs.", call. = FALSE)
  }
  taxon_names <- unique(trimws(as.character(taxon_names)))
  taxon_names <- taxon_names[!is.na(taxon_names) & nzchar(taxon_names)]
  progress <- NULL
  on.exit(if (!is.null(progress)) close(progress), add = TRUE)
  has_api_key <- !is.null(api_key) && nzchar(api_key)
  request_interval <- if (has_api_key) 0.11 else 0.34
  last_request <- -Inf
  request_ncbi <- function(url, body, stage) {
    for (attempt in seq_len(10L)) {
      wait <- last_request + request_interval - as.numeric(Sys.time())
      if (wait > 0) Sys.sleep(wait)
      last_request <<- as.numeric(Sys.time())
      handle <- curl::new_handle()
      curl::handle_setopt(handle, post = TRUE, postfields = body)
      response <- tryCatch(curl::curl_fetch_memory(url, handle = handle),
                           error = function(e) e)
      html_error <- !inherits(response, "error") &&
        grepl("^\\s*<(?:html|!doctype html)",
              rawToChar(head(response$content, 100L)),
              ignore.case = TRUE, perl = TRUE)
      if (!inherits(response, "error") && response$status_code == 200L &&
          !html_error) {
        return(response)
      }
      if (!inherits(response, "error") &&
          !html_error &&
          !response$status_code %in% c(429L, 500L, 502L, 503L, 504L)) {
        break
      }
    }
    detail <- if (inherits(response, "error")) conditionMessage(response) else
      paste0("HTTP ", response$status_code, ": ",
             substr(rawToChar(response$content), 1L, 300L))
    stop(stage, " failed after ", attempt, " attempt(s): ", detail,
         call. = FALSE)
  }
  terms <- character()
  organism_terms <- c(
    paste0("txid", taxids, "[Organism:exp]"),
    sprintf('"%s"[Organism:exp]', taxon_names)
  )
  if (length(organism_terms)) {
    terms <- c(terms, paste0("(", paste(organism_terms, collapse = " OR "), ")"))
  }
  if (!is.null(gene)) terms <- c(terms, sprintf("(%s[Gene] OR %s[Title] OR %s[All Fields])", gene, gene, gene))
  if (!is.null(date_from) || !is.null(date_to)) {
    if (is.null(date_from)) date_from <- "1900/01/01"
    if (is.null(date_to)) date_to <- format(Sys.Date(), "%Y/%m/%d")
    terms <- c(terms, sprintf("%s:%s[PDAT]", date_from, date_to))
  }
  if (!is.null(min_length) || !is.null(max_length)) {
    search_min_length <- if (is.null(min_length)) 1 else as.integer(min_length)
    search_max_length <- if (is.null(max_length)) 1000000000 else as.integer(max_length)
    terms <- c(
      terms,
      sprintf("%d:%d[SLEN]", search_min_length, search_max_length)
    )
  }
  term <- paste(c(terms, ncbi_query), collapse = " AND ")
  species_count <- if (!is.null(taxon_count)) {
    as.integer(taxon_count)
  } else if (length(taxids) || length(taxon_names)) {
    length(taxids) + length(taxon_names)
  } else {
    NA_integer_
  }
  if (is.na(species_count)) {
    show_time("Searching NCBI Entrez...")
  } else {
    show_time(paste0("Searching NCBI Entrez for ", species_count,
                     " species..."))
  }

  search_params <- c(
    db = "nuccore", term = term, retmode = "json",
    retstart = "0", retmax = "10000",
    tool = "METAeDNA"
  )
  if (!is.null(email)) search_params <- c(search_params, email = email)
  if (has_api_key) search_params <- c(search_params, api_key = api_key)
  search_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  search_body <- paste(
    names(search_params),
    vapply(search_params, utils::URLencode, "", reserved = TRUE),
    sep = "=", collapse = "&"
  )
  search_response <- request_ncbi(search_url, search_body, "NCBI ESearch")
  search <- jsonlite::fromJSON(rawToChar(search_response$content))$esearchresult
  record_count <- suppressWarnings(as.numeric(search$count))
  if (length(record_count) != 1L || is.na(record_count)) {
    stop("NCBI ESearch did not return a record count: ",
         substr(rawToChar(search_response$content), 1L, 500L), call. = FALSE)
  }
  records_to_fetch <- if (is.infinite(max_records)) {
    record_count
  } else {
    min(record_count, as.numeric(max_records))
  }
  if (!is.finite(records_to_fetch) || records_to_fetch <= 0) {
    warning("No NCBI records matched the query.", call. = FALSE)
    return(Biostrings::DNAStringSet())
  }

  if (length(batch_size) != 1L || is.na(batch_size) ||
      !is.numeric(batch_size) || batch_size < 1) {
    stop("batch_size must be a positive number.", call. = FALSE)
  }
  summary_batch_size <- 10000L
  fetch_batch_size <- as.integer(min(batch_size, 500L))
  expressions <- as.character(
    regex_expression[!is.na(regex_expression) & nzchar(regex_expression)]
  )
  show_time(paste0("Screening ", records_to_fetch, " NCBI record titles..."))
  path <- if (is.null(output_fasta)) {
    tempfile(fileext = ".fasta")
  } else {
    dir.create(dirname(output_fasta), recursive = TRUE, showWarnings = FALSE)
    tempfile(pattern = paste0(basename(output_fasta), ".part-"),
             tmpdir = dirname(output_fasta))
  }
  if (is.null(output_fasta)) on.exit(unlink(path), add = TRUE)
  taxid_map <- list()
  missing_summaries <- 0L
  position <- 0L
  progress <- utils::txtProgressBar(min = 0, max = records_to_fetch, width = 40, style = 3)
  while (position < records_to_fetch) {
    i <- length(taxid_map) + 1L
    page_size <- min(summary_batch_size, records_to_fetch - position)
    if (position == 0L) {
      ids <- as.character(search$idlist)
    } else {
      search_params["retstart"] <- as.character(position)
      search_params["retmax"] <- as.character(page_size)
      search_body <- paste(
        names(search_params),
        vapply(search_params, utils::URLencode, "", reserved = TRUE),
        sep = "=", collapse = "&"
      )
      page_response <- request_ncbi(search_url, search_body,
                                    paste0("NCBI ESearch batch ", i))
      page <- jsonlite::fromJSON(rawToChar(page_response$content))$esearchresult
      ids <- as.character(page$idlist)
    }
    ids <- head(ids, page_size)
    if (length(ids) != page_size) {
      stop("NCBI ESearch batch ", i, " returned ", length(ids),
           " of ", page_size, " requested UIDs. Response: ",
           substr(rawToChar(if (position == 0L) search_response$content
                            else page_response$content), 1L, 500L),
           call. = FALSE)
    }
    summary_params <- c(
      db = "nuccore", id = paste(ids, collapse = ","),
      retmax = as.character(page_size),
      retmode = "xml", version = "2.0", tool = "METAeDNA"
    )
    if (!is.null(email)) summary_params <- c(summary_params, email = email)
    if (has_api_key) summary_params <- c(summary_params, api_key = api_key)
    summary_body <- paste(
      names(summary_params),
      vapply(summary_params, utils::URLencode, "", reserved = TRUE),
      sep = "=", collapse = "&"
    )
    summary_response <- request_ncbi(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
      summary_body, paste0("NCBI ESummary batch ", i)
    )
    summary_xml <- tryCatch(
      xml2::read_xml(summary_response$content, options = "NONET"),
      error = function(e) NULL
    )
    summary_set <- if (is.null(summary_xml)) NULL else
      xml2::xml_find_first(summary_xml, ".//DocumentSummarySet")
    summary_status <- if (is.null(summary_set)) NA_character_ else
      xml2::xml_attr(summary_set, "status")
    records <- if (is.null(summary_set)) NULL else
      xml2::xml_find_all(summary_set, "./DocumentSummary")
    if (!identical(summary_status, "OK") || !length(records)) {
      stop("NCBI ESummary batch ", i, " did not return a valid summary set: ",
           gsub("[[:space:]]+", " ",
                substr(rawToChar(summary_response$content), 1L, 500L)),
           call. = FALSE)
    }
    missing_summaries <- missing_summaries +
      length(setdiff(ids, xml2::xml_attr(records, "uid")))
    taxid_map[[i]] <- data.frame(
      uid = xml2::xml_attr(records, "uid"),
      accession = xml2::xml_text(xml2::xml_find_first(records, "./AccessionVersion")),
      taxid = xml2::xml_text(xml2::xml_find_first(records, "./TaxId")),
      title = xml2::xml_text(xml2::xml_find_first(records, "./Title"))
    )
    taxid_map[[i]] <- taxid_map[[i]][
      !is.na(taxid_map[[i]]$accession) & nzchar(taxid_map[[i]]$accession) &
        !is.na(taxid_map[[i]]$taxid) & nzchar(taxid_map[[i]]$taxid),
      , drop = FALSE
    ]
    if (length(expressions) && !is.null(taxid_map[[i]])) {
      keep <- rep(FALSE, nrow(taxid_map[[i]]))
      for (expression in expressions) {
        keep <- keep | grepl(expression, taxid_map[[i]]$title,
                             ignore.case = TRUE, perl = TRUE)
      }
      taxid_map[[i]] <- taxid_map[[i]][keep, , drop = FALSE]
    }
    position <- position + page_size
    utils::setTxtProgressBar(progress, position)
  }
  close(progress)
  progress <- NULL
  if (missing_summaries) {
    warning("NCBI did not return summaries for ", missing_summaries,
            " searched UIDs; those records were skipped.", call. = FALSE)
  }
  taxid_map <- do.call(rbind, taxid_map)
  if (is.null(taxid_map) || !nrow(taxid_map)) {
    stop("No NCBI records matched regex_expression and the length filters.",
         call. = FALSE)
  }
  taxid_map <- taxid_map[!duplicated(taxid_map$accession), , drop = FALSE]
  fetch_starts <- seq.int(1L, nrow(taxid_map), by = fetch_batch_size)
  show_time(paste0("Matched ", nrow(taxid_map), " accessions; downloading FASTA in ",
                   length(fetch_starts), " batches..."))
  progress <- utils::txtProgressBar(min = 0, max = nrow(taxid_map), width = 40, style = 3)
  for (i in seq_along(fetch_starts)) {
    selected_accessions <- taxid_map$accession[seq.int(
      fetch_starts[i], min(fetch_starts[i] + fetch_batch_size - 1L,
                           nrow(taxid_map))
    )]
    fetch_params <- c(
      db = "nuccore", id = paste(selected_accessions, collapse = ","),
      rettype = "fasta", retmode = "text", tool = "METAeDNA"
    )
    if (!is.null(email)) fetch_params <- c(fetch_params, email = email)
    if (has_api_key) fetch_params <- c(fetch_params, api_key = api_key)
    fetch_body <- paste(
      names(fetch_params),
      vapply(fetch_params, utils::URLencode, "", reserved = TRUE),
      sep = "=", collapse = "&"
    )
    fetch_response <- request_ncbi(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
      fetch_body, paste0("NCBI EFetch batch ", i)
    )
    fasta_text <- rawToChar(fetch_response$content)
    if (!startsWith(fasta_text, ">")) {
      stop("NCBI EFetch batch ", i, " did not return FASTA: ",
           substr(fasta_text, 1L, 300L), call. = FALSE)
    }
    cat(fasta_text, "\n", file = path, append = TRUE, sep = "")
    utils::setTxtProgressBar(progress, min(fetch_starts[i] + length(selected_accessions) - 1L,
                                           nrow(taxid_map)))
  }
  close(progress)
  progress <- NULL

  seqs <- Biostrings::readDNAStringSet(path)
  seqs <- clean_gb_reference_sequences(seqs, min_length, max_length)
  if (!length(seqs)) {
    stop("No downloaded sequences matched the length and title filters.",
         call. = FALSE)
  }
  taxid_map <- taxid_map[
    !is.na(taxid_map$accession) & grepl("^[1-9][0-9]*$", taxid_map$taxid) &
      taxid_map$accession %in% parse_gb_accession(names(seqs)),
    , drop = FALSE
  ]
  if (!nrow(taxid_map)) {
    stop("No downloaded accession has an NCBI TaxID.", call. = FALSE)
  }
  attr(seqs, "taxid_map") <- unique(taxid_map[, c("accession", "taxid")])
  if (!is.null(output_fasta)) {
    Biostrings::writeXStringSet(seqs, output_fasta)
    unlink(path)
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
#' @param n_threads BLAST threads.
#' @param min_length,max_length Optional sequence length filters.
#' @param sql_path Optional local NCBI taxonomy SQLite database.
#' @param taxon_names Optional taxon whitelist applied through \code{sql_path}.
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
#'   n_threads = 8,
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
                             n_threads = 1, min_length = NULL, max_length = NULL,
                             sql_path = NULL, taxon_names = NULL, taxidlist = NULL) {
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
        "-evalue", evalue, "-num_threads", n_threads,
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

    if (!is.null(sql_path) && !is.null(taxon_names) && length(round_sequences)) {
      ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
      taxid <- get_taxid(sql_path, parse_gb_accession(names(round_sequences)))
      taxonomy <- get_taxonomy(sql_path, unique(stats::na.omit(taxid$taxid)),
                               rank_name = ranks)
      taxonomy <- dplyr::left_join(taxid, taxonomy, by = dplyr::join_by(taxid == query_taxid))
      keep <- apply(taxonomy[, intersect(ranks, names(taxonomy)), drop = FALSE], 1, function(x) any(taxon_names %in% x))
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

#' Build GenBank reference sequences online, locally, or from seeds
#'
#' Online mode searches NCBI Entrez, checks record titles, and downloads matching
#' FASTA sequences. Local mode selects records from a prepared nt BLAST database
#' by TaxID and title. Mine mode expands seed sequences by iterative BLAST
#' against a local database. Set `makeblastdb = TRUE` to format the resulting
#' FASTA as a BLAST database.
#'
#' @param mode One of \code{"online"}, \code{"local"}, or \code{"mine"}.
#' @param taxon_names NCBI taxon names. In local extraction mode, `ncbi_db` is
#'   required to resolve names to TaxIDs.
#' @param taxids NCBI TaxIDs or a table with an `ncbi_taxid` column.
#' @param gene,ncbi_query,date_from,date_to Entrez search settings for online
#'   mode. `ncbi_query` uses NCBI syntax; `regex_expression` is an optional
#'   additional title filter.
#' @param regex_expression Title regular expression(s) used by both extraction
#'   modes. In online mode, titles are checked before FASTA download.
#' @param seed_fasta Seed FASTA required for mine mode.
#' @param database Local BLAST database prefix for local and mine modes.
#' @param output_dir Output directory.
#' @param email,api_key NCBI contact/API settings. Online requests are limited
#'   to 3 per second without an API key or 10 per second with one.
#' @param max_records,batch_size Online retrieval controls.
#' @param max_rounds,max_target_seqs,round_seed_limit,perc_identity,qcov_hsp_perc,evalue
#'   Iterative BLAST controls for mine mode.
#' @param n_threads Number of local BLAST threads. Online requests are sequential.
#' @param min_length,max_length Optional sequence length filters.
#' @param ncbi_db Local NCBI taxonomy SQLite database for resolving names in
#'   local extraction mode.
#' @param sql_path Optional taxonomy SQLite database for mine-mode validation.
#' @param blast_exe_path Optional BLAST+ executable directory.
#' @param makeblastdb Build a BLAST database from the resulting FASTA.
#' @param database_name Basename of the FASTA and BLAST database.
#'
#' @return A list containing `sequences`, `fasta`, `blast_database`, and
#'   source-specific results. `taxid_map` is the accession-to-TaxID file,
#'   written whenever TaxIDs are available, including when `makeblastdb = FALSE`.
#'   It is `NULL` when no mapping is available.
#' @examples
#' \dontrun{
#' build_gb_reference(
#'   mode = "local", database = "blast_database/nt", taxids = 41944,
#'   regex_expression = marker_regex$trnL,
#'   output_dir = "reference/local_trnL", makeblastdb = TRUE
#' )
#' build_gb_reference(
#'   mode = "online", taxids = 41944,
#'   ncbi_query = marker_entrez_query$trnL,
#'   regex_expression = marker_regex$trnL,
#'   output_dir = "reference/online_trnL", makeblastdb = TRUE
#' )
#' build_gb_reference(
#'   mode = "mine", seed_fasta = "reference/trnL_seeds.fasta",
#'   database = "blast_database/nt", output_dir = "reference/mined_trnL",
#'   max_rounds = 3, makeblastdb = TRUE
#' )
#' }
#' @export
build_gb_reference <- function(mode = c("online", "local", "mine"),
                       taxon_names = NULL, taxids = NULL,
                       gene = NULL, ncbi_query = NULL,
                       regex_expression = NULL,
                       date_from = NULL, date_to = NULL,
                       seed_fasta = NULL, database = NULL,
                       output_dir = "build_gb_reference",
                       email = NULL, api_key = NULL,
                       max_records = 500, batch_size = 100,
                       max_rounds = 3, max_target_seqs = 50,
                       round_seed_limit = 200,
                       perc_identity = 85, qcov_hsp_perc = 70,
                       evalue = "1e-20", n_threads = 1,
                       min_length = NULL, max_length = NULL,
                       ncbi_db = NULL, sql_path = NULL,
                       blast_exe_path = NULL, makeblastdb = FALSE,
                       database_name = "reference") {
  mode <- match.arg(mode)
  if (mode == "mine" && is.null(seed_fasta)) {
    stop("seed_fasta is required for mine mode.", call. = FALSE)
  }
  if (mode != "mine" && !is.null(seed_fasta)) {
    stop("seed_fasta is only used in mine mode.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  start_time <- Sys.time()

  if (mode == "online") {
    fasta <- file.path(output_dir, paste0(database_name, ".fasta"))
    seqs <- fetch_gb_reference_online(
      taxon_names = taxon_names, taxids = taxids,
      gene = gene, ncbi_query = ncbi_query,
      regex_expression = regex_expression,
      date_from = date_from, date_to = date_to,
      max_records = max_records, batch_size = batch_size,
      email = email, api_key = api_key, n_threads = n_threads,
      min_length = min_length, max_length = max_length,
      output_fasta = fasta
    )
    taxid_pairs <- attr(seqs, "taxid_map")
    source_result <- seqs
  } else if (mode == "local") {
    source_result <- fetch_gb_reference_local(
      database = database, output_dir = output_dir,
      taxids = taxids, taxon_names = taxon_names, ncbi_db = ncbi_db,
      regex_expression = regex_expression,
      min_length = min_length, max_length = max_length,
      output_fasta = paste0(database_name, ".fasta"),
      blast_exe_path = blast_exe_path
    )
    fasta <- source_result$fasta
    if (is.na(fasta)) stop("No local nt records matched the filters.", call. = FALSE)
    seqs <- Biostrings::readDNAStringSet(fasta)
    taxid_pairs <- source_result$records[, c("accession", "taxid"), drop = FALSE]
  } else {
    source_result <- mine_gb_reference_local(
      seed_fasta = seed_fasta, database = database, output_dir = output_dir,
      blast_exe_path = blast_exe_path, max_rounds = max_rounds,
      max_target_seqs = max_target_seqs, round_seed_limit = round_seed_limit,
      perc_identity = perc_identity, qcov_hsp_perc = qcov_hsp_perc,
      evalue = evalue, n_threads = n_threads,
      min_length = min_length, max_length = max_length,
      sql_path = sql_path, taxon_names = taxon_names
    )
    fasta <- source_result$output_fasta
    seqs <- source_result$sequences
    taxid_pairs <- NULL
  }

  taxid_map <- NULL
  blast_database <- NULL
  if (!is.null(taxid_pairs) && nrow(taxid_pairs)) {
    taxid_pairs <- taxid_pairs[
      !is.na(taxid_pairs$taxid) &
        grepl("^[1-9][0-9]*$", as.character(taxid_pairs$taxid)),
      , drop = FALSE
    ]
    if (nrow(taxid_pairs)) {
      taxid_map <- file.path(output_dir, "taxid_map.tsv")
      utils::write.table(unique(taxid_pairs), taxid_map,
                         sep = "\t", quote = FALSE,
                         row.names = FALSE, col.names = FALSE)
    }
  }
  if (isTRUE(makeblastdb)) {
    blast_database <- file.path(output_dir, database_name)
    status <- blast_cmd(
      base_cmd = "makeblastdb",
      args = c("-in", fasta, "-dbtype", "nucl", "-parse_seqids",
               if (!is.null(taxid_map)) c("-taxid_map", taxid_map),
               "-out", blast_database),
      blast_exe_path = blast_exe_path
    )
    if (!identical(as.integer(status), 0L)) {
      stop("makeblastdb failed for the reference FASTA.", call. = FALSE)
    }
  }
  result <- list(mode = mode, sequences = seqs, fasta = fasta,
                 output_fasta = fasta, blast_database = blast_database,
                 taxid_map = taxid_map, source = source_result)
  if (mode == "mine") {
    result$hits <- source_result$hits
    result$rounds <- source_result$rounds
  }
  show_time(secs = as.numeric(Sys.time() - start_time, units = "secs"), mode = 1)
  result
}

#' Extract annotated nucleotide records from a local nt database
#'
#' Use NCBI `blastdbcmd` to list accessions belonging to selected TaxIDs, filter
#' them by title and sequence length, and then extract only the selected records
#' as FASTA. The BLAST database directory is added temporarily to `BLASTDB` so
#' that `taxdb.btd` and `taxdb.bti` can be found during taxonomic filtering.
#'
#' @param database Prefix of a local NCBI nucleotide BLAST database, such as
#'   `"blastdb/nt"`.
#' @param output_dir Directory used for the TaxID list, record table,
#'   accession list, and output FASTA.
#' @param taxids NCBI TaxIDs or a table containing `ncbi_taxid`, such as the
#'   result of `col_ncbi_crosswalk()`.
#' @param crosswalk Optional result returned by `col_ncbi_crosswalk()`.
#'   Existing calls remain supported; new calls can pass this table as `taxids`.
#' @param taxon_names Optional NCBI taxon name or vector of names. Names are matched
#'   exactly against the local taxonomy database, including its synonyms when
#'   `name_aliases` is available. Descendant sequences are included.
#' @param ncbi_db Local NCBI taxonomy SQLite database required for `taxon_names`.
#' @param regex_expression Optional text or regular expressions matched against
#'   nt record titles without regard to letter case, for example
#'   `c("\\brbcL\\b", "ribulose-1,5-bisphosphate carboxylase")`. A record is
#'   retained when it matches any expression. When omitted, every record for
#'   the selected TaxIDs is extracted.
#' @param min_length,max_length Optional sequence-length limits.
#' @param output_fasta Name of the resulting FASTA file inside `output_dir`.
#' @param blast_exe_path Optional BLAST+ executable directory.
#'
#' @return A list containing the selected record table and paths to the TaxID
#'   list, accession list, and FASTA file. The record table contains accession,
#'   sequence length, and title columns.
#'
#' @details The standard nt BLAST database stores sequence titles but not the
#' full structured GenBank feature table. Consequently, `regex_expression`
#' filtering is based on the title reported by `blastdbcmd`. Ordinary text such
#' as `"rbcL"` performs a substring match. Expressions containing regular
#' expression syntax, such as `"\\brbcL\\b"` or `"rbcL|matK"`, use that syntax
#' during matching.
#'
#' @examples
#' \dontrun{
#' names <- col_ncbi_crosswalk(
#'   c("Quercus robur", "Quercus ilex"),
#'   col_db = "col.sqlite",
#'   ncbi_db = "ncbi_taxonomy.sqlite"
#' )
#'
#' local_reference <- fetch_gb_reference_local(
#'   database = "blastdb/nt",
#'   output_dir = "quercus_rbcl",
#'   taxids = names,
#'   regex_expression = marker_regex$rbcL,
#'   min_length = 400,
#'   max_length = 2000
#' )
#' local_reference$fasta
#'
#' fetch_gb_reference_local(
#'   database = "blastdb/nt", output_dir = "quercus_rbcl",
#'   taxon_names = c("Quercus robur", "Quercus ilex"),
#'   ncbi_db = "ncbi_taxonomy.sqlite",
#'   regex_expression = marker_regex$rbcL
#' )
#' }
#' @export
fetch_gb_reference_local <- function(database,
                                     output_dir,
                                     taxids = NULL,
                                     crosswalk = NULL,
                                     regex_expression = NULL,
                                     min_length = NULL,
                                     max_length = NULL,
                                     output_fasta = "reference.fasta",
                                     blast_exe_path = NULL,
                                     taxon_names = NULL,
                                     ncbi_db = NULL) {
  .check_path(database, "database")
  .check_path(output_dir, "output_dir")
  .check_path(output_fasta, "output_fasta")
  if (!is.null(min_length)) .check_number(min_length, "min_length", 0, .Machine$integer.max, TRUE)
  if (!is.null(max_length)) .check_number(max_length, "max_length", 0, .Machine$integer.max, TRUE)
  if (is.data.frame(taxids) && !"ncbi_taxid" %in% names(taxids)) {
    stop("taxids table must contain an ncbi_taxid column.", call. = FALSE)
  }
  if (!is.null(crosswalk) && !"ncbi_taxid" %in% names(crosswalk)) {
    stop("crosswalk must contain an ncbi_taxid column.", call. = FALSE)
  }
  taxids <- c(
    if (is.data.frame(taxids)) suppressWarnings(as.integer(taxids$ncbi_taxid))
    else suppressWarnings(as.integer(taxids)),
    if (!is.null(crosswalk)) suppressWarnings(as.integer(crosswalk$ncbi_taxid))
  )
  if (length(taxon_names)) {
    if (is.null(ncbi_db)) {
      stop("ncbi_db is required when taxon names are supplied.", call. = FALSE)
    }
    con <- DBI::dbConnect(RSQLite::SQLite(), ncbi_db)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    name_table <- if ("name_aliases" %in% DBI::dbListTables(con)) "name_aliases" else "names"
    taxon_names <- unique(trimws(as.character(taxon_names)))
    taxon_names <- taxon_names[!is.na(taxon_names) & nzchar(taxon_names)]
    matched <- DBI::dbGetQuery(con, sprintf(
      "SELECT DISTINCT tax_id, name FROM %s WHERE name IN (%s)",
      DBI::dbQuoteIdentifier(con, name_table),
      paste(DBI::dbQuoteString(con, taxon_names), collapse = ",")
    ))
    missing_names <- setdiff(taxon_names, matched$name)
    if (length(missing_names)) {
      warning("NCBI taxonomy did not match: ",
              paste(missing_names, collapse = ", "), call. = FALSE)
    }
    taxids <- c(taxids, matched$tax_id)
  }
  taxids <- unique(taxids)
  taxids <- taxids[!is.na(taxids) & taxids > 0L]
  if (!length(taxids)) {
    stop("Provide a valid TaxID, taxon name, or table with NCBI matches.",
         call. = FALSE)
  }
  show_time(paste0("Searching local nt for ", length(taxids), " taxa..."))

  if (!is.null(min_length) && (!is.numeric(min_length) || min_length < 0)) {
    stop("min_length must be a non-negative number.", call. = FALSE)
  }
  if (!is.null(max_length) && (!is.numeric(max_length) || max_length < 0)) {
    stop("max_length must be a non-negative number.", call. = FALSE)
  }
  if (!is.null(min_length) && !is.null(max_length) && min_length > max_length) {
    stop("min_length cannot be greater than max_length.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  taxid_file <- file.path(output_dir, "taxids.txt")
  record_file <- file.path(output_dir, "nt_records.tsv")
  accession_file <- file.path(output_dir, "selected_accessions.txt")
  fasta_file <- file.path(output_dir, output_fasta)
  writeLines(as.character(taxids), taxid_file)

  database_dir <- normalizePath(dirname(database), mustWork = FALSE)
  old_blastdb <- Sys.getenv("BLASTDB", unset = NA_character_)
  blastdb_paths <- c(database_dir, old_blastdb[!is.na(old_blastdb) & nzchar(old_blastdb)])
  Sys.setenv(BLASTDB = paste(unique(blastdb_paths), collapse = .Platform$path.sep))
  on.exit({
    if (is.na(old_blastdb)) {
      Sys.unsetenv("BLASTDB")
    } else {
      Sys.setenv(BLASTDB = old_blastdb)
    }
  }, add = TRUE)

  status <- blast_cmd(
    base_cmd = "blastdbcmd",
    args = c(
      "-db", shQuote(database),
      "-taxidlist", shQuote(taxid_file),
      "-target_only",
      "-outfmt", shQuote("%a\t%l\t%S\t%T\t%t"),
      "-out", shQuote(record_file)
    ),
    blast_exe_path = blast_exe_path
  )
  if (!identical(as.integer(status), 0L) || !file.exists(record_file)) {
    stop("blastdbcmd could not list records for the supplied TaxIDs.",
         call. = FALSE)
  }

  if (file.info(record_file)$size == 0) {
    records <- data.frame(
      accession = character(), length = numeric(),
      scientific_name_GenBank = character(), taxid = numeric(),
      title = character()
    )
  } else {
    records <- utils::read.delim(
      record_file,
      header = FALSE,
      sep = "\t",
      quote = "",
      col.names = c("accession", "length", "scientific_name_GenBank","taxid","title"),
      colClasses = c("character", "numeric", "character", "numeric" ,"character"),
      fill = TRUE,
      comment.char = "",
      stringsAsFactors = FALSE
    )
  }
  keep <- rep(TRUE, nrow(records))
  expressions <- as.character(
    regex_expression[!is.na(regex_expression) & nzchar(regex_expression)]
  )
  if (length(expressions)) {
    keep_regex <- rep(FALSE, nrow(records))
    for (expression in expressions) {
      keep_regex <- keep_regex |
        grepl(expression, records$title, ignore.case = TRUE, perl = TRUE)
    }
    keep <- keep & keep_regex
  }
  if (!is.null(min_length)) keep <- keep & records$length >= min_length
  if (!is.null(max_length)) keep <- keep & records$length <= max_length

  selected <- which(keep & !duplicated(records$accession))
  records <- records[selected, , drop = FALSE]
  utils::write.table(records, file.path(output_dir, "selected_records.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  if (!nrow(records)) {
    show_time("Selected 0 sequences from 0 taxa.")
    writeLines(character(), accession_file)
    return(list(
      records = records,
      taxid_file = taxid_file,
      accession_file = accession_file,
      fasta = NA_character_
    ))
  }

  writeLines(records$accession, accession_file)
  status <- blast_cmd(
    base_cmd = "blastdbcmd",
    args = c(
      "-db", shQuote(database),
      "-entry_batch", shQuote(accession_file),
      "-target_only",
      "-outfmt", shQuote("%f"),
      "-out", shQuote(fasta_file)
    ),
    blast_exe_path = blast_exe_path
  )
  if (!identical(as.integer(status), 0L) || !file.exists(fasta_file)) {
    stop("blastdbcmd could not extract the selected sequences.", call. = FALSE)
  }
  show_time(paste0("Extracted ", nrow(records), " sequences from ",
                   length(unique(records$taxid)), " taxa."))
  unlink(record_file)
  list(
    records = records,
    taxid_file = taxid_file,
    accession_file = accession_file,
    fasta = fasta_file
  )
}
