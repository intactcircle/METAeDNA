#' Run the METAeDNA browser interface
#'
#' Open a local Shiny interface for FASTQ quality inspection, sequence processing,
#' DADA2, blank correction, replicate merging, classification and community matrices.
#' Each step runs in a separate R process and saves its parameters, log and results.
#'
#' @param work_dir Directory for analysis results. Relative input paths are
#'   resolved against this directory. The default is the current directory.
#' @param launch.browser Open the application in a browser.
#' @param host Address on which to listen. The default allows local access only.
#' @param port Optional TCP port; Shiny selects one when NULL.
#'
#' @details Install the optional packages `shiny` and `callr` before use.
#'   File paths refer to the machine running R. Choose existing FASTQ, FASTA,
#'   CSV/TSV or RDS files; large sequencing files do not need browser uploads.
#'   Run steps individually and select their saved files as subsequent inputs.
#'   Quality previews inspect the beginning of each file. Overlap mismatch
#'   rates are available for merged FASTQ from `pairend()` and describe retained
#'   overlaps, not the error rate of all raw reads.
#'   The sequence viewer reads sequences from the current result or a local
#'   FASTA file. It displays an AliView-style, base-coloured sequence panel;
#'   click a row to open the complete sequence and adjust its font size.
#'   The interface permits local file access and is intended for a trusted
#'   workstation, not an unauthenticated public server.
#' @return The value returned by `shiny::runApp()` when the application stops.
#' @importFrom stats setNames
#' @importFrom utils head
#' @examples
#' \dontrun{
#' install.packages(c("shiny", "callr"))
#' run_shiny()
#' }
#' @export
run_shiny <- function(work_dir = ".", launch.browser = interactive(),
                      host = "127.0.0.1", port = NULL) {
  for (package in c("shiny", "callr")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Install package '", package, "' to use run_shiny().", call. = FALSE)
  }
  .check_path(work_dir, "work_dir")
  .check_flag(launch.browser, "launch.browser")
  .check_path(host, "host")
  if (!is.null(port)) .check_number(port, "port", 1, 65535, TRUE)
  if (!dir.exists(work_dir) && !dir.create(work_dir, recursive = TRUE)) stop("Cannot create work_dir.", call. = FALSE)
  shiny::runApp(.metaedna_app(normalizePath(work_dir)), launch.browser = launch.browser, host = host, port = port)
}

# Field defaults also define the input type; no R expressions are evaluated.
.workflow_fields <- function() {
  list(
    quality = list(fastq_1 = "", fastq_2 = "", max_records = 100000),
    pairend = list(fastq_1 = "", fastq_2 = "", min_overlap = 15, max_mismatches = 3, n_threads = 1, batch_size = 100000),
    reorder_fastq = list(fastq_1 = "", fastq_2 = ""),
    demultiplex = list(fastq_file = "", reaction_table = "", primer_mismatch = 3, with_tag = TRUE, n_threads = 1),
    barcode_split = list(fastq_1 = "", fastq_2 = "", reaction_table = "", with_primer = TRUE, n_threads = 1, batch_size = 100000),
    dada2 = list(barcode_files = "", experiment = "run1", trunc_len = c(0, 0), max_ee = c(2, 2), n_threads = 1),
    collapse_sequences = list(fasta_file = "", n_threads = 1, batch_size = 100000),
    denoise = list(fasta_file = "", sequence_min_L = 1, sequence_max_L = 1000, min_count = 1,
                   threshold = 0.05, max_diff = 1, max_indel = 0, keep_internal = FALSE, n_threads = 1),
    vsearch_remove_chimeras = list(fasta_file = "", mode = c("uchime_denovo", "uchime2_denovo", "uchime3_denovo", "uchime_ref"),
                                  ref_fasta = "", vsearch_exe_path = "", n_threads = 1),
    vsearch_cluster = list(fasta_file = "", mode = c("cluster_fast", "cluster_size", "cluster_smallmem", "cluster_unoise"),
                           identity = 0.97, minsize = 8, unoise_alpha = 2, n_threads = 1, vsearch_exe_path = ""),
    fasta_to_table = list(fasta_file = "", n_threads = 1),
    decontamination = list(sequence_table = "", control_table = "", method = c("maximum", "mean", "threshold", "both_threshold"),
                           control_names = "FB,EB,PB", threshold = 10, missing_columns = c("zero", "error")),
    merge_replicate = list(sequence_table = "", control_table = "", method = c("mean_rra", "max_rra", "both_occur"),
                           repeat_times = 1, replicate_number = "", present_absent = FALSE),
    blastn_cmd = list(fasta_file = "", database = "", blast_exe_path = "", perc_identity = 90, qcov_hsp_perc = 90, max_target_seqs = 50, n_threads = 1),
    get_taxonomy = list(sql_path = "", blast_result = "", col_path = ""),
    find_lca = list(df = "", sequences = "", identity_threshold = ""),
    merge_to_matrix = list(merged_table = "", control_table = "", taxonomy_table = "", taxa_in_rows = FALSE),
    rarefaction_extrapolation = list(sequence_table = "", taxonomy_table = "", bin_size = "", endpoint = "", n_threads = 1,
                                    show_observed = TRUE, show_extrapolation = TRUE),
    in_silico = list(fasta_file = "", forward = "", reverse = "", mis_match = 3, try_reverse = TRUE, trim = TRUE),
    makeblastdb = list(fasta_file = "", blast_exe_path = ""),
    download_blast_db = list(database_path = "blast_download", resume = TRUE, md5_threads = 1),
    prepare_blast_database = list(file_path = "blast_download", n_threads = 1),
    build_ncbi_taxonomy_database = list(taxonomy = TRUE, accession = FALSE, create_index = TRUE, chunk_size = 1000000),
    build_col_taxonomy_database = list(col_data = ""),
    fetch_gb_reference_local = list(database = "", taxids = "", taxon_names = "", ncbi_db = "", regex_expression = "", min_length = "", max_length = "", blast_exe_path = "")
  )
}

.workflow_step_help <- c(
  quality = "Preview Phred quality and base composition before filtering.",
  pairend = "Merge paired reads from R1 and R2 using their best accepted overlap.",
  reorder_fastq = "Reorder R2 so that its read identifiers follow R1.",
  demultiplex = "Assign merged reads to PCR reactions using primers and reaction tags.",
  barcode_split = "Split tagged paired-end reads into one forward/reverse file pair per PCR reaction.",
  dada2 = "Filter, denoise, merge and remove chimeras with DADA2.",
  collapse_sequences = "Collapse identical sequences and retain per-reaction counts.",
  denoise = "Label or remove low-abundance neighbours of abundant sequences.",
  vsearch_remove_chimeras = "Remove chimeric sequences with VSEARCH.",
  vsearch_cluster = "Cluster sequences with a selected VSEARCH algorithm.",
  fasta_to_table = "Convert processed FASTA metadata into a sequence count table.",
  decontamination = "Correct sample counts using the assigned blank reactions.",
  merge_replicate = "Normalize and merge PCR replicates within each experiment.",
  blastn_cmd = "Search query sequences against a nucleotide BLAST database.",
  get_taxonomy = "Attach NCBI taxonomy and optionally reconcile names with COL.",
  find_lca = "Calculate the lowest common ancestor from a BLAST taxonomy table.",
  merge_to_matrix = "Aggregate classified sequences into a community matrix.",
  rarefaction_extrapolation = "Calculate reaction-level rarefaction and extrapolation curves.",
  in_silico = "Extract amplicons matching a primer pair from reference sequences.",
  makeblastdb = "Build a local nucleotide BLAST database from FASTA.",
  download_blast_db = "Download and verify the NCBI nt database archives.",
  prepare_blast_database = "Extract downloaded BLAST archives into a usable database.",
  build_ncbi_taxonomy_database = "Build the local NCBI taxonomy SQLite database.",
  build_col_taxonomy_database = "Build a local Catalogue of Life SQLite database.",
  fetch_gb_reference_local = "Retrieve taxon- and annotation-matched records from local nt."
)

.workflow_parameter_help <- c(
  fastq_1 = "Forward (R1) FASTQ path. Plain text and gzip files are accepted.",
  fastq_2 = "Reverse (R2) FASTQ path. Leave empty only when the selected step permits it.",
  fastq_file = "Input FASTQ path.", fasta_file = "Input FASTA path.",
  reaction_table = "CSV or TSV file containing reaction tags and primers.",
  control_table = "CSV, TSV or RDS reaction design table.", sequence_table = "CSV, TSV or RDS sequence count table.",
  taxonomy_table = "CSV, TSV or RDS taxonomy table joined through SEQ.", merged_table = "Merged sequence table path.",
  n_threads = "Number of worker threads or processes; use 1 for low-memory execution.",
  batch_size = "Records read per batch. Larger batches use more memory.", max_records = "Maximum FASTQ records inspected; 0 reads the complete file.",
  min_overlap = "Minimum accepted overlap between R1 and reverse-complemented R2.",
  max_mismatches = "Maximum mismatches allowed in the selected overlap.", primer_mismatch = "Maximum primer mismatches; tag matching remains exact.",
  with_tag = "Require tag matching as well as primer matching.", with_primer = "Remove matched primers after splitting.",
  trunc_len = "R1,R2 truncation lengths, for example 220,180; 0 keeps the full read.", max_ee = "Maximum expected errors for R1,R2, for example 2,2.",
  threshold = "Cutoff used by the selected correction or denoising method.", identity_threshold = "Optional rank thresholds, for example species=100,genus=98,family=95.",
  control_names = "Comma-separated blank categories, for example FB,EB,PB.", missing_columns = "Fill missing reaction columns with zero or stop with an error.",
  replicate_number = "Optional fixed averaging denominator; empty uses detected replicate count.", endpoint = "Optional extrapolation end depth; empty uses 110% of the maximum depth.",
  sequence_min_L = "Minimum retained sequence length.", sequence_max_L = "Maximum retained sequence length.", min_count = "Minimum sequence abundance.",
  max_diff = "Maximum substitutions for denoising neighbours.", max_indel = "Maximum short indels for denoising neighbours.",
  keep_internal = "Keep sequences labelled internal; otherwise their counts are removed.", database = "BLAST database prefix, not an individual archive.",
  taxids = "One or more comma-separated NCBI TaxIDs.", taxon_names = "Comma-separated NCBI taxon names; requires ncbi_db.", regex_expression = "Plain annotation text or a regular expression used to select accessions.",
  blast_exe_path = "Optional BLAST+ executable directory; empty searches PATH.", vsearch_exe_path = "Optional VSEARCH executable path; empty searches PATH."
)

.workflow_required <- list(
  quality = "fastq_1", pairend = c("fastq_1", "fastq_2"), reorder_fastq = c("fastq_1", "fastq_2"),
  demultiplex = c("fastq_file", "reaction_table"), barcode_split = c("fastq_1", "fastq_2", "reaction_table"),
  dada2 = "barcode_files", collapse_sequences = "fasta_file", denoise = "fasta_file",
  vsearch_remove_chimeras = "fasta_file", vsearch_cluster = "fasta_file", fasta_to_table = "fasta_file",
  decontamination = c("sequence_table", "control_table"), merge_replicate = c("sequence_table", "control_table"),
  blastn_cmd = c("fasta_file", "database"), get_taxonomy = c("sql_path", "blast_result"), find_lca = "df",
  merge_to_matrix = c("merged_table", "control_table", "taxonomy_table"), rarefaction_extrapolation = "sequence_table",
  in_silico = c("fasta_file", "forward", "reverse"), makeblastdb = "fasta_file",
  download_blast_db = "database_path", prepare_blast_database = "file_path",
  build_col_taxonomy_database = "col_data", fetch_gb_reference_local = "database"
)

.read_workflow_table <- function(path) {
  .check_path(path, "Table path")
  if (tolower(tools::file_ext(path)) == "rds") return(readRDS(path))
  separator <- if (grepl("\t", readLines(path, n = 1, warn = FALSE), fixed = TRUE)) "\t" else ","
  utils::read.table(path, header = TRUE, sep = separator, quote = "\"", comment.char = "",
                    check.names = FALSE, stringsAsFactors = FALSE)
}

.run_workflow_step <- function(step, arguments, output_dir) {
  saveRDS(list(step = step, arguments = arguments, started = Sys.time()), file.path(output_dir, "parameters.rds"))
  start <- Sys.time()
  paths <- character()
  result <- NULL
  if (step == "quality") {
    files <- unlist(arguments[c("fastq_1", "fastq_2")])
    files <- files[nzchar(files)]
    if (!length(files)) stop("Provide at least one FASTQ path.")
    result <- lapply(files, fastq_profile, max_records = arguments$max_records)
    for (i in seq_along(result)) {
      result[[i]]$positions$read <- names(files)[i]
      result[[i]]$summary$read <- names(files)[i]
    }
    result <- list(positions = do.call(rbind, lapply(result, `[[`, "positions")),
                   summary = do.call(rbind, lapply(result, `[[`, "summary")))
    utils::write.csv(result$positions, file.path(output_dir, "quality_by_position.csv"), row.names = FALSE)
    utils::write.csv(result$summary, file.path(output_dir, "quality_summary.csv"), row.names = FALSE)
  } else if (step == "dada2") {
    files <- .read_workflow_table(arguments$barcode_files)
    .check_columns(files, c("sample", "forward", "reverse"), "barcode_files")
    result <- dada2_pipeline(files$forward, files$reverse, file.path(output_dir, "dada2"),
                             sample_names = files$sample, trunc_len = arguments$trunc_len,
                             max_ee = arguments$max_ee, n_threads = arguments$n_threads)
    table <- as.data.frame(t(result$sequence_table), check.names = FALSE)
    names(table) <- paste0("count_", rownames(result$sequence_table))
    table$SEQ <- colnames(result$sequence_table)
    table$ID <- table$SEQ_ID <- sprintf("ASV_%07d", seq_len(nrow(table)))
    table$experiment <- arguments$experiment
    table$total_count <- rowSums(t(result$sequence_table))
    saveRDS(table, file.path(output_dir, "sequence_table.rds"))
    paths <- file.path(output_dir, "sequence_table.rds")
  } else if (step == "makeblastdb") {
    db <- file.path(output_dir, "reference")
    status <- blast_cmd("makeblastdb", c("-in", shQuote(arguments$fasta_file), "-dbtype", "nucl", "-parse_seqids", "-out", shQuote(db)), arguments$blast_exe_path)
    if (status != 0) stop("makeblastdb failed: ", status)
    result <- list(database = db)
  } else {
    for (name in intersect(names(arguments), c("sequence_table", "control_table", "taxonomy_table", "merged_table", "df"))) {
      if (!is.null(arguments[[name]])) arguments[[name]] <- .read_workflow_table(arguments[[name]])
    }
    if (step == "find_lca") {
      sequences <- arguments$sequences
      arguments$sequences <- NULL
    }
    file_outputs <- switch(step,
      pairend = c(out_file = "merged.fastq.gz"),
      reorder_fastq = c(output_fastq = "reordered_R2.fastq.gz"),
      demultiplex = c(identified_output = "identified.fasta.gz", unidentified_output = "unidentified.fasta.gz"),
      collapse_sequences = c(output_fasta = "collapsed.fasta.gz"),
      denoise = c(output_fasta = "denoised.fasta.gz"),
      vsearch_remove_chimeras = c(output_fasta = "nonchimera.fasta"),
      blastn_cmd = c(output = "blast.tsv"), character())
    for (name in names(file_outputs)) arguments[[name]] <- file.path(output_dir, file_outputs[[name]])
    if (step %in% c("barcode_split", "fetch_gb_reference_local")) arguments$output_dir <- file.path(output_dir, "files")
    if (step == "vsearch_cluster") arguments$output_path <- file.path(output_dir, "clusters")
    if (step == "prepare_blast_database") arguments$database_path <- file.path(output_dir, "database")
    if (step == "build_ncbi_taxonomy_database") {
      arguments$taxdump_dir <- file.path(output_dir, "taxdump")
      arguments$sql_path <- file.path(output_dir, "ncbi.sqlite")
    }
    if (step == "build_col_taxonomy_database") arguments$db_path <- file.path(output_dir, "col.sqlite")
    result <- do.call(get(step, envir = asNamespace("METAeDNA")), arguments)
    paths <- unname(unlist(arguments[names(file_outputs)]))
    if (length(paths) && any(!file.exists(paths))) stop("A required output file was not created.")
    if (step == "pairend") {
      result <- fastq_profile(paths[[1]], max_records = 0)
      result$positions$read <- "merged"
      result$summary$read <- "merged"
    }
    if (step == "find_lca" && nzchar(sequences)) {
      sequence_table <- .read_workflow_table(sequences)
      .check_columns(sequence_table, c("SEQ_ID", "SEQ"), "sequences")
      result$SEQ <- sequence_table$SEQ[match(result$qseqid, sequence_table$SEQ_ID)]
    }
    if (step == "in_silico") {
      path <- file.path(output_dir, "amplicons.fasta")
      Biostrings::writeXStringSet(result, path)
      paths <- c(paths, path)
    }
    if (step == "barcode_split") paths <- c(paths, file.path(arguments$output_dir, "barcode_files.csv"))
  }
  result_path <- file.path(output_dir, "result.rds")
  saveRDS(result, result_path)
  if (is.data.frame(result) || is.matrix(result)) {
    utils::write.csv(result, file.path(output_dir, "result.csv"), row.names = is.matrix(result))
  }
  if (is.data.frame(result) && all(c("SEQ_ID", "SEQ") %in% names(result))) {
    sequences <- unique(result[c("SEQ_ID", "SEQ")])
    dna <- Biostrings::DNAStringSet(sequences$SEQ)
    names(dna) <- sequences$SEQ_ID
    Biostrings::writeXStringSet(dna, file.path(output_dir, "queries.fasta"))
  }
  list(result_path = result_path, files = unique(c(paths, list.files(output_dir, full.names = TRUE))),
       elapsed = as.numeric(difftime(Sys.time(), start, units = "secs")))
}

.metaedna_app <- function(work_dir) {
  steps <- .workflow_fields()
  ui <- shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(paste0(
      "body{background:#f7f8fa}.well{background:white}pre{max-height:320px;overflow:auto}",
      ".container-fluid{max-width:1500px}.btn-primary{background:#3B4992}h2{color:#24304a}",
      ".sequence-view{overflow:auto;max-height:650px;background:#fff;border:1px solid #ddd;font-family:monospace;white-space:nowrap}",
      ".sequence-row{display:flex;min-width:max-content;border-bottom:1px solid #f0f0f0;cursor:pointer;line-height:1.55}",
      ".sequence-row:hover{background:#eef2ff}.sequence-name{position:sticky;left:0;z-index:1;width:180px;overflow:hidden;",
      "text-overflow:ellipsis;background:inherit;border-right:1px solid #ddd;padding:0 8px;color:#222}",
      ".sequence-bases{padding-left:8px;letter-spacing:1px}.base-A{color:#087f23}.base-C{color:#1565c0}",
      ".base-G{color:#e17b00}.base-T{color:#c62828}.base-N{color:#777}.base-gap{color:#aaa}",
      ".sequence-detail{font-family:monospace;overflow-wrap:anywhere;line-height:1.7;letter-spacing:1px}",
      ".help-label{cursor:help;border-bottom:1px dotted #777}.step-help{min-height:42px;color:#4d596b}",
      ".run-status{padding:8px 10px;margin-bottom:8px;background:#fff;border-left:4px solid #3B4992}"
    ))),
    shiny::titlePanel("METAeDNA - Sequencing workflow"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(width = 4, shiny::helpText(paste("Working directory:", work_dir)),
        shiny::tags$div(title = "Choose one processing stage. Hover over its parameters for details.",
                        shiny::selectInput("step", "Analysis step", choices = names(steps))),
        shiny::uiOutput("step_help"),
        shiny::actionButton("example", "Use packaged FASTQ example"),
        shiny::helpText("Paths refer to this computer. Table inputs accept CSV, TSV or RDS. Run one step, then use its output in the next step."),
        shiny::uiOutput("fields"), shiny::actionButton("run", "Run step", class = "btn-primary"),
        shiny::actionButton("cancel", "Cancel"), shiny::hr(),
        shiny::selectInput("view_job", "View a completed run", choices = character()),
        shiny::uiOutput("artifact_controls"),
        shiny::downloadButton("download_result", "Download result RDS")),
      shiny::mainPanel(width = 8,
        shiny::tags$div(class = "run-status", shiny::textOutput("status")),
        shiny::tabsetPanel(
          shiny::tabPanel("Quality", shiny::helpText("Phred+33: mean quality with the 10-90% range. Previews use the first records, not a random sample."), shiny::plotOutput("quality"), shiny::tableOutput("qc_summary")),
          shiny::tabPanel("Base composition", shiny::plotOutput("bases")),
          shiny::tabPanel("Sequence viewer",
            shiny::fluidRow(shiny::column(9, shiny::textInput("sequence_file", "FASTA file", placeholder = "Optional local FASTA path")),
                            shiny::column(3, shiny::br(), shiny::actionButton("load_sequences", "Load FASTA"))),
            shiny::helpText("Sequences from the current result are shown automatically. Click a row to inspect the complete sequence."),
            shiny::uiOutput("sequence_summary"), shiny::uiOutput("sequence_viewer")),
          shiny::tabPanel("Results", shiny::tableOutput("preview"), shiny::plotOutput("result_plot")),
          shiny::tabPanel("Files and history", shiny::tableOutput("history"), shiny::tableOutput("files")),
          shiny::tabPanel("Log", shiny::verbatimTextOutput("log")))))
  )
  server <- function(input, output, session) {
    run_files <- file.path(list.dirs(work_dir, recursive = FALSE, full.names = TRUE), "run.rds")
    run_files <- run_files[file.exists(run_files)]
    previous_jobs <- lapply(run_files, function(path) tryCatch(readRDS(path), error = function(e) NULL))
    previous_jobs <- Filter(function(job) is.list(job) && !is.null(job$step) && !is.null(job$status), previous_jobs)
    previous_jobs <- lapply(previous_jobs, function(job) {
      if (identical(job$status, "running")) job$status <- "interrupted"
      job
    })
    state <- shiny::reactiveValues(process = NULL, jobs = previous_jobs, selected = NULL,
                                   sequences = NULL, log = "", status = "Ready")

    sequence_data <- shiny::reactive({
      if (!is.null(state$sequences)) return(state$sequences)
      value <- state$selected
      if (inherits(value, "DNAStringSet")) {
        sequence <- as.character(value)
        return(data.frame(name = if (is.null(names(value))) paste0("sequence_", seq_along(sequence)) else names(value),
                          sequence = sequence, stringsAsFactors = FALSE))
      }
      if (is.data.frame(value) && "SEQ" %in% names(value)) {
        name_column <- intersect(c("SEQ_ID", "ID", "qseqid", "scientific_name"), names(value))
        name <- if (length(name_column)) as.character(value[[name_column[1]]]) else paste0("sequence_", seq_len(nrow(value)))
        sequences <- unique(data.frame(name = name, sequence = as.character(value$SEQ), stringsAsFactors = FALSE))
        return(sequences[!is.na(sequences$sequence) & nzchar(sequences$sequence), , drop = FALSE])
      }
      if (is.list(value) && is.matrix(value$sequence_table)) {
        sequence <- colnames(value$sequence_table)
        return(data.frame(name = paste0("ASV_", seq_along(sequence)), sequence = sequence, stringsAsFactors = FALSE))
      }
      NULL
    })

    base_spans <- function(sequence, limit = Inf) {
      sequence <- toupper(sequence)
      if (is.finite(limit)) sequence <- substr(sequence, 1, limit)
      bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
      classes <- ifelse(bases %in% c("A", "C", "G", "T", "N"), paste0("base-", bases),
                        ifelse(bases == "-", "base-gap", "base-N"))
      shiny::tagList(Map(function(base, class) shiny::tags$span(class = class, base), bases, classes))
    }

    shiny::observeEvent(input$load_sequences, {
      tryCatch({
        .check_path(input$sequence_file, "FASTA file")
        if (!file.exists(input$sequence_file)) stop("FASTA file not found.")
        dna <- Biostrings::readDNAStringSet(input$sequence_file)
        sequence <- as.character(dna)
        state$sequences <- data.frame(
          name = if (is.null(names(dna))) paste0("sequence_", seq_along(sequence)) else names(dna),
          sequence = sequence, stringsAsFactors = FALSE
        )
      }, error = function(e) shiny::showNotification(conditionMessage(e), type = "error"))
    })
    output$fields <- shiny::renderUI({
      fields <- steps[[input$step]]
      shiny::tagList(lapply(names(fields), function(name) {
        value <- fields[[name]]
        id <- paste0("arg_", name)
        help <- unname(.workflow_parameter_help[name])
        if (is.na(help)) help <- paste("Value passed to", paste0(input$step, "()."))
        label <- shiny::tags$span(class = "help-label", title = help, name)
        control <- if (is.logical(value)) shiny::checkboxInput(id, label, value)
        else if (is.character(value) && length(value) > 1) shiny::selectInput(id, label, value)
        else if (is.numeric(value) && length(value) == 1) shiny::numericInput(id, label, value)
        else shiny::textInput(id, label, paste(value, collapse = ","))
        shiny::tags$div(title = help, control)
      }))
    })
    output$step_help <- shiny::renderUI({
      shiny::tags$div(class = "step-help", title = .workflow_step_help[[input$step]],
                      shiny::tags$strong(input$step), shiny::tags$br(), .workflow_step_help[[input$step]])
    })
    output$artifact_controls <- shiny::renderUI({
      files <- unique(unlist(lapply(state$jobs, function(job) {
        if (!identical(job$status, "complete") || is.null(job$result)) return(character())
        job$result$files
      })))
      files <- files[file.exists(files)]
      fields <- steps[[input$step]]
      targets <- names(fields)[vapply(fields, function(value) is.character(value) && length(value) == 1L, logical(1))]
      if (!length(files) || !length(targets)) return(NULL)
      shiny::tagList(
        shiny::selectInput("artifact_file", shiny::tags$span(class = "help-label", title = "Output file from a completed step.", "Previous output"), files),
        shiny::selectInput("artifact_target", shiny::tags$span(class = "help-label", title = "Parameter that will receive the selected file path.", "Use as parameter"), targets),
        shiny::actionButton("use_artifact", "Use selected file", title = "Copy this path into the selected parameter.")
      )
    })
    shiny::observeEvent(input$use_artifact, {
      shiny::req(nzchar(input$artifact_file), nzchar(input$artifact_target))
      shiny::updateTextInput(session, paste0("arg_", input$artifact_target), value = input$artifact_file)
    })
    shiny::observeEvent(input$example, {
      values <- c(fastq_1 = system.file("extdata", "trnL_40k_1.fq.gz", package = "METAeDNA"),
                  fastq_2 = system.file("extdata", "trnL_40k_2.fq.gz", package = "METAeDNA"),
                  reaction_table = system.file("extdata", "trnL_table.txt", package = "METAeDNA"),
                  control_table = system.file("extdata", "control_table.csv", package = "METAeDNA"))
      for (name in intersect(names(values), names(steps[[input$step]]))) {
        shiny::updateTextInput(session, paste0("arg_", name), value = values[[name]])
      }
    })
    shiny::observeEvent(input$run, {
      if (!is.null(state$process) && state$process$is_alive()) {
        shiny::showNotification("A step is already running.", type = "warning")
        return()
      }
      tryCatch({
        step <- input$step
        fields <- steps[[step]]
        arguments <- setNames(lapply(names(fields), function(name) input[[paste0("arg_", name)]]), names(fields))
        for (name in names(arguments)) {
          value <- arguments[[name]]
          if (is.null(value)) stop("Wait for the parameter fields to load.")
          if (name %in% c("trunc_len", "max_ee") || (name == "taxids" && nzchar(value))) {
            value <- suppressWarnings(as.numeric(trimws(strsplit(value, ",", fixed = TRUE)[[1]])))
            if (!length(value) || any(!is.finite(value))) stop(name, " must contain only comma-separated numbers.")
          }
          if (name == "taxon_names" && nzchar(value)) {
            value <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
          }
          if (name %in% c("taxids", "taxon_names", "ncbi_db") && identical(value, "")) value <- NULL
          if (name == "control_names") {
            value <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
            if (!length(value) || any(!nzchar(value))) stop("control_names contains an empty category.")
          }
          if (name == "identity_threshold") {
            if (!nzchar(value)) value <- NULL else {
              parts <- strsplit(strsplit(value, ",", fixed = TRUE)[[1]], "=", fixed = TRUE)
              if (any(lengths(parts) != 2L)) stop("Use species=100,genus=98,family=95 for identity_threshold.")
              value <- setNames(suppressWarnings(as.numeric(vapply(parts, `[`, "", 2))), trimws(vapply(parts, `[`, "", 1)))
              if (any(!is.finite(value)) || any(!nzchar(names(value)))) stop("identity_threshold contains an invalid rank or value.")
            }
          }
          if (name %in% c("replicate_number", "bin_size", "endpoint", "min_length", "max_length")) {
            value <- if (nzchar(value)) suppressWarnings(as.numeric(value)) else NULL
            if (length(value) && !is.finite(value)) stop(name, " must be a number or left empty.")
          }
          if (name %in% c("ref_fasta", "vsearch_exe_path", "blast_exe_path", "col_path", "regex_expression", "taxonomy_table") && identical(value, "")) value <- NULL
          arguments[name] <- list(value)
        }
        required <- .workflow_required[[step]]
        absent <- required[vapply(arguments[required], function(value) {
          is.null(value) || !length(value) || (is.character(value) && all(!nzchar(trimws(value))))
        }, logical(1))]
        if (length(absent)) stop("Required field(s): ", paste(absent, collapse = ", "), ".")
        directory <- tempfile(paste0(step, "_"), tmpdir = work_dir)
        if (!dir.create(directory)) stop("Cannot create run directory.")
        state$log <- file.path(directory, "run.log")
        manifest <- list(step = step, status = "running", directory = directory, started = Sys.time())
        state$jobs <- c(state$jobs, list(manifest))
        state$status <- paste("Running", step)
        saveRDS(manifest, file.path(directory, "run.rds"))
        package_path <- normalizePath(find.package("METAeDNA"), mustWork = TRUE)
        source_mode <- dir.exists(file.path(package_path, "src"))
        if (source_mode && !requireNamespace("pkgload", quietly = TRUE)) {
          stop("Install pkgload to run background steps from a development source tree.")
        }
        library_paths <- if (source_mode) .libPaths() else unique(c(dirname(package_path), .libPaths()))
        state$process <- callr::r_bg(function(step, arguments, directory, work_dir,
                                             package_path, source_mode) {
          setwd(work_dir)
          if (source_mode) {
            pkgload::load_all(package_path, quiet = TRUE)
          } else {
            loadNamespace("METAeDNA")
          }
          utils::getFromNamespace(".run_workflow_step", "METAeDNA")(step, arguments, directory)
        }, args = list(step, arguments, directory, work_dir, package_path, source_mode),
        libpath = library_paths,
        stdout = state$log, stderr = state$log, supervise = TRUE)
      }, error = function(e) {
        jobs <- state$jobs
        if (length(jobs) && identical(jobs[[length(jobs)]]$status, "running")) {
          jobs[[length(jobs)]]$status <- "failed"
          jobs[[length(jobs)]]$error <- conditionMessage(e)
          jobs[[length(jobs)]]$finished <- Sys.time()
          saveRDS(jobs[[length(jobs)]], file.path(jobs[[length(jobs)]]$directory, "run.rds"))
          state$jobs <- jobs
        }
        state$status <- conditionMessage(e)
        shiny::showNotification(state$status, type = "error", duration = NULL)
      })
    })
    shiny::observe({
      shiny::invalidateLater(1000, session)
      process <- state$process
      if (is.null(process) || process$is_alive()) return()
      i <- length(state$jobs)
      jobs <- state$jobs
      jobs[[i]]$finished <- Sys.time()
      result <- tryCatch(process$get_result(), error = identity)
      if (inherits(result, "error")) {
        log_lines <- if (file.exists(state$log)) readLines(state$log, warn = FALSE) else character()
        log_lines <- gsub("\033\\[[0-9;]*m", "", log_lines, perl = TRUE)
        useful <- trimws(log_lines)
        useful <- useful[nzchar(useful) & useful != "Execution halted" &
                           !grepl("^[-]+$|^Error in doTryCatch", useful)]
        detail <- if (length(useful)) tail(useful, 1L) else conditionMessage(result)
        jobs[[i]]$status <- "failed"
        jobs[[i]]$error <- detail
        state$status <- paste("Failed:", detail)
        shiny::showNotification(state$status, type = "error", duration = NULL)
      } else {
        jobs[[i]]$status <- "complete"
        jobs[[i]]$result <- result
        state$selected <- readRDS(result$result_path)
        state$sequences <- NULL
        state$status <- sprintf("Completed %s in %.1f seconds", jobs[[i]]$step, result$elapsed)
      }
      saveRDS(jobs[[i]], file.path(jobs[[i]]$directory, "run.rds"))
      state$jobs <- jobs
      state$process <- NULL
    })
    shiny::observe({
      jobs <- state$jobs
      choices <- which(vapply(jobs, function(job) identical(job$status, "complete") &&
                                !is.null(job$result) && file.exists(job$result$result_path), logical(1)))
      labels <- if (length(choices)) vapply(jobs[choices], function(job) {
        paste(job$step, basename(job$directory), sep = " - ")
      }, character(1)) else character()
      shiny::updateSelectInput(session, "view_job", choices = setNames(as.character(choices), labels),
                               selected = if (length(choices)) tail(as.character(choices), 1L) else character())
    })
    shiny::observeEvent(input$cancel, {
      if (!is.null(state$process) && state$process$is_alive()) {
        state$process$kill_tree()
        jobs <- state$jobs
        i <- length(jobs)
        jobs[[i]]$status <- "cancelled"
        jobs[[i]]$finished <- Sys.time()
        saveRDS(jobs[[i]], file.path(jobs[[i]]$directory, "run.rds"))
        state$jobs <- jobs
        state$process <- NULL
        state$status <- "Cancelled; partial outputs remain in the run directory."
      }
    })
    session$onSessionEnded(function() {
      process <- shiny::isolate(state$process)
      if (!is.null(process) && process$is_alive()) process$kill_tree()
    })
    shiny::observeEvent(input$view_job, {
      i <- suppressWarnings(as.integer(input$view_job))
      if (length(i) == 1 && !is.na(i) && i <= length(state$jobs)) {
        job <- state$jobs[[i]]
        if (!is.null(job$result) && file.exists(job$result$result_path)) {
          state$selected <- readRDS(job$result$result_path)
          state$sequences <- NULL
        }
      }
    })
    output$status <- shiny::renderText(state$status)
    output$log <- shiny::renderText({
      shiny::invalidateLater(1000, session)
      if (nzchar(state$log) && file.exists(state$log)) paste(tail(readLines(state$log, warn = FALSE), 100), collapse = "\n") else "No log yet."
    })
    output$history <- shiny::renderTable({
      if (!length(state$jobs)) return(NULL)
      do.call(rbind, lapply(state$jobs, function(x) data.frame(step = x$step, status = x$status, directory = x$directory)))
    })
    output$files <- shiny::renderTable({
      if (!length(state$jobs)) return(NULL)
      data.frame(file = unlist(lapply(state$jobs, function(x) list.files(x$directory, recursive = TRUE, full.names = TRUE))))
    })
    output$qc_summary <- shiny::renderTable({
      value <- state$selected
      if (is.list(value) && !is.null(value$positions)) value$summary
    }, digits = 5)
    output$quality <- shiny::renderPlot({
      value <- state$selected
      shiny::req(is.list(value), !is.null(value$positions), nrow(value$positions) > 0)
      d <- value$positions
      ggplot2::ggplot(d, ggplot2::aes(x = .data$position, y = .data$mean_quality, colour = .data$read, fill = .data$read)) +
        ggplot2::geom_boxplot(ggplot2::aes(group = .data$position, ymin = .data$q10,
          lower = .data$q25, middle = .data$median, upper = .data$q75, ymax = .data$q90),
          stat = "identity", width = 0.7, alpha = 0.3, linewidth = 0.3) +
        ggplot2::geom_line() + ggplot2::facet_wrap(~read) +
        ggplot2::labs(x = "Base position", y = "Phred quality") +
        ggplot2::scale_colour_manual(values = c("#3B4992", "#EE0000")) +
        ggplot2::scale_fill_manual(values = c("#3B4992", "#EE0000")) + ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(strip.background = ggplot2::element_blank())
    })
    output$bases <- shiny::renderPlot({
      value <- state$selected
      shiny::req(is.list(value), !is.null(value$positions), nrow(value$positions) > 0)
      d <- tidyr::pivot_longer(value$positions, cols = dplyr::all_of(c("A", "C", "G", "T", "N_other")), names_to = "base", values_to = "frequency")
      ggplot2::ggplot(d, ggplot2::aes(.data$position, .data$frequency, colour = .data$base)) + ggplot2::geom_line() +
        ggplot2::facet_wrap(~read) + ggplot2::scale_colour_manual(values = c("#3B4992", "#EE0000", "#008B45", "#631879", "#008280")) +
        ggplot2::labs(x = "Base position", y = "Base frequency") + ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(strip.background = ggplot2::element_blank())
    })
    output$sequence_summary <- shiny::renderUI({
      sequences <- sequence_data()
      if (is.null(sequences) || !nrow(sequences)) return(shiny::helpText("No sequences are available in the current result."))
      lengths <- nchar(sequences$sequence)
      shiny::helpText(sprintf("%s sequences; length %s-%s bp. The alignment view shows at most 200 sequences and 2,000 bases per row.",
                              format(nrow(sequences), big.mark = ","), min(lengths), max(lengths)))
    })
    output$sequence_viewer <- shiny::renderUI({
      sequences <- sequence_data()
      shiny::req(!is.null(sequences), nrow(sequences) > 0)
      shown <- seq_len(min(200L, nrow(sequences)))
      rows <- lapply(shown, function(i) {
        shiny::tags$div(
          class = "sequence-row",
          onclick = sprintf("Shiny.setInputValue('selected_sequence', %d, {priority: 'event'})", i),
          shiny::tags$span(class = "sequence-name", title = sequences$name[i], sequences$name[i]),
          shiny::tags$span(class = "sequence-bases", base_spans(sequences$sequence[i], 2000L))
        )
      })
      shiny::tags$div(class = "sequence-view", rows)
    })
    shiny::observeEvent(input$selected_sequence, {
      sequences <- sequence_data()
      index <- as.integer(input$selected_sequence)
      shiny::req(!is.null(sequences), !is.na(index), index >= 1L, index <= nrow(sequences))
      shiny::showModal(shiny::modalDialog(
        title = sequences$name[index], size = "l", easyClose = TRUE,
        shiny::helpText(paste("Length:", nchar(sequences$sequence[index]), "bp")),
        shiny::sliderInput("sequence_font_size", "Font size", min = 10, max = 28, value = 16, step = 1),
        shiny::uiOutput("sequence_detail"), footer = shiny::modalButton("Close")
      ))
    })
    output$sequence_detail <- shiny::renderUI({
      sequences <- sequence_data()
      index <- as.integer(input$selected_sequence)
      shiny::req(!is.null(sequences), !is.na(index), index >= 1L, index <= nrow(sequences))
      font_size <- if (is.null(input$sequence_font_size)) 16 else input$sequence_font_size
      shiny::tags$div(class = "sequence-detail", style = paste0("font-size:", font_size, "px"),
                      base_spans(sequences$sequence[index]))
    })
    shiny::outputOptions(output, "sequence_detail", suspendWhenHidden = FALSE)
    output$preview <- shiny::renderTable({
      value <- state$selected
      if (is.data.frame(value) || is.matrix(value)) head(value[, seq_len(min(12, ncol(value))), drop = FALSE], 20)
      else if (is.list(value) && !is.null(value$tracking)) value$tracking
    })
    output$result_plot <- shiny::renderPlot({
      value <- state$selected
      if (is.list(value) && !is.null(value$fig)) print(value$fig)
    })
    output$download_result <- shiny::downloadHandler(filename = function() "METAeDNA_result.rds",
      content = function(file) { shiny::req(!is.null(state$selected)); saveRDS(state$selected, file) })
  }
  shiny::shinyApp(ui, server)
}
