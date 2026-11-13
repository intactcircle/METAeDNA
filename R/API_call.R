# Vsearch -------------------------------------
## vsearch command --------------------------------
#' Execute a VSEARCH command with automatic executable detection
#'
#' \code{vsearch_cmd()} locates VSEARCH and runs the supplied arguments using
#' \code{system2()}. The executable path and command are printed before execution.
#'
#' @param vsearch_exe_path (character or NULL) optional path to the VSEARCH
#'   executable. If \code{NULL}, the function searches for \code{"vsearch"} in
#'   the system PATH using \code{detect_cmd()}. If a directory is provided, the
#'   function attempts to locate the binary inside that directory.
#'
#' @param args (character vector) command-line arguments passed directly to
#'   \code{system2()}. Each element is treated as one argument; for example:
#'   \code{c("--threads", "8", "--cluster_fast",
#'   "sequences.fasta")}.
#'
#' @return
#' Invisibly returns the exit status of \code{system2()}, typically \code{0} for
#' successful execution. The function also prints:
#' \itemize{
#'   \item the resolved path to the VSEARCH executable,
#'   \item the complete assembled command string before execution.
#' }
#'
#' @importFrom stringr str_detect
#'
#' @examples
#' \dontrun{
#' vsearch_cmd(args = "--version")
#' }#'
#' @export
vsearch_cmd <-function(vsearch_exe_path = NULL, args){
  # detect path -------------------
  # show_time(paste0("Searching path..."))
  exe <- detect_cmd("vsearch", path = vsearch_exe_path)
  show_time(paste0("Executable file found in \"", exe, "\""))
  # make command -----------------------
  show_time(paste("COMMAND:", exe, paste(args, collapse = " "), "\n"))
  res <- system2(command = exe, args = args)
  invisible(res)
}

## vsearch command --------------------------------
#' Remove chimeric sequences using VSEARCH uchime-denovo or uchime-ref algorithms
#'
#' \code{vsearch_remove_chimeras()} runs VSEARCH's \code{uchime_denovo},
#' \code{uchime2_denovo}, \code{uchime3_denovo}, or \code{uchime_ref} algorithm.
#' Reference-based detection requires a reference FASTA file.
#'
#' @param vsearch_exe_path (character or NULL) optional path to the VSEARCH
#'   executable. If \code{NULL}, the function searches for \code{"vsearch"} on
#'   the system PATH using \code{detect_cmd()}.
#'
#' @param fasta_file (character) path to the input FASTA file containing the
#'   sequences to be checked for chimeras.
#'
#' @param output_fasta (character) path to the output FASTA file. The content
#'   depends on \code{output_mode} and may contain non-chimeras, chimeras, or
#'   uchime alignment outputs.
#'
#' @param mode (character) chimera-detection algorithm to use. One of:
#'   \itemize{
#'     \item \code{"uchime_denovo"}
#'     \item \code{"uchime2_denovo"}
#'     \item \code{"uchime3_denovo"}
#'     \item \code{"uchime_ref"} (requires \code{ref_fasta})
#'   }
#'   Passed to VSEARCH as \code{--uchime_denovo}, \code{--uchime2_denovo}, etc.
#'
#' @param output_mode (character) output format produced by VSEARCH. One of:
#'   \itemize{
#'     \item \code{"nonchimeras"} — write only non-chimeric sequences
#'     \item \code{"chimeras"} — write detected chimeras
#'     \item \code{"uchimealns"} — write UCHIME alignment information
#'     \item \code{"uchimeout"} — write detailed chimera scoring output
#'   }
#'   Passed to VSEARCH as \code{--nonchimeras}, \code{--chimeras}, etc.
#'
#' @param ref_fasta (character or NULL) reference FASTA file required only when
#'   \code{mode = "uchime_ref"}. For all reference-free modes, this argument may
#'   be left \code{NULL}.
#'
#' @param minuniquesize (integer or NULL) optional minimum abundance threshold
#'   (VSEARCH \code{--minuniquesize}) used mainly in de novo chimera detection.
#'   If \code{NULL}, the option is omitted.
#'
#' @param n_threads (integer or NULL) number of CPU threads to pass to VSEARCH via
#'   \code{--threads}. If \code{NULL}, threading is not explicitly set.
#'
#' @return
#' Returns the exit status of the executed VSEARCH command (integer). Also prints:
#' \itemize{
#'   \item the resolved VSEARCH executable path,
#'   \item the full constructed command string,
#'   \item runtime messages from \code{vsearch_cmd()}.
#' }
#'
#' @examples
#' \dontrun{
#' # collapsed.fasta is the output of collapse_sequences().
#' vsearch_remove_chimeras(
#'   fasta_file = "collapsed.fasta", output_fasta = "nonchimera.fasta",
#'   mode = "uchime3_denovo", n_threads = 4
#' )
#' }#'
#' @export
vsearch_remove_chimeras <- function(
    vsearch_exe_path = NULL,
    fasta_file,
    output_fasta,
    mode = c("uchime_denovo", "uchime2_denovo", "uchime3_denovo", "uchime_ref"),
    output_mode = c("nonchimeras", "chimeras" ,   "uchimealns", "uchimeout"),
    ref_fasta = NULL,
    minuniquesize = NULL,
    n_threads = NULL
){
  mode <- match.arg(mode)
  output_mode <- match.arg(output_mode)
  .check_path(fasta_file, "fasta_file")
  .check_path(output_fasta, "output_fasta")
  if (!is.null(n_threads)) .check_number(n_threads, "n_threads", 1, .Machine$integer.max, TRUE)
  if (!is.null(minuniquesize)) .check_number(minuniquesize, "minuniquesize", 1, .Machine$integer.max, TRUE)
  if((mode ==  "uchime_ref") & is.null(ref_fasta)){stop("In uchime_ref mode, a specific reference database (ref_fasta) must be provided as input.")}
  mode <- paste0("--", mode)
  output_mode <- paste0("--", output_mode)
  if(is.null(ref_fasta)){
    args <- c(mode, shQuote(fasta_file), output_mode, shQuote(output_fasta),
              if (!is.null(n_threads)) c("--threads", n_threads))
  }else{
    args <- c(mode, shQuote(fasta_file), output_mode, shQuote(output_fasta), "--db", shQuote(ref_fasta),
              if (!is.null(n_threads)) c("--threads", n_threads))
  }
  if (!is.null(minuniquesize)) args <- c(args, "--minuniquesize", minuniquesize)
  status <- vsearch_cmd(vsearch_exe_path = vsearch_exe_path, args = args)
  if (status != 0L) stop("VSEARCH failed with exit status ", status, ".", call. = FALSE)
  invisible(status)
}

## vsearch_cluster ----------------------------------
#' Cluster a FASTA file or sequence table using VSEARCH
#'
#' The function accepts either a FASTA file through \code{fasta_file} or a
#' METAeDNA sequence-abundance table through \code{sequence_table}. Both input modes
#' write cluster centroids and a \code{.uc} mapping file. When \code{sequence_table}
#' is supplied, cluster membership is additionally merged back into the table
#' and numeric abundance columns are summed for each experiment and centroid.
#'
#' @param vsearch_exe_path (character or NULL) optional path to the VSEARCH
#'   executable. If \code{NULL}, the function searches for \code{"vsearch"}
#'   via \code{detect_cmd()}.
#'
#' @param fasta_file (character) path to the FASTA file containing the input
#'   sequences. Supply either \code{fasta_file} or \code{sequence_table}.
#' @param sequence_table (data.frame or tibble) optional sequence-abundance table
#'   containing \code{SEQ_ID}, \code{SEQ}, and \code{total_count}. Supply either
#'   \code{sequence_table} or \code{fasta_file}.
#'
#' @param output_path (character) directory where clustering output files will
#'   be written. The function creates:
#'   \itemize{
#'     \item \code{vsearch_cluster.fas} — cluster centroid sequences,
#'     \item \code{vsearch_cluster.uc}  — VSEARCH cluster mapping file.
#'   }
#' @param identity (numeric) sequence identity threshold passed to VSEARCH via
#'   \code{--id}. Typical values range from 0.90 to 0.99, depending on the desired
#'   clustering level (e.g., OTU97 uses \code{0.97}).
#'
#' @param mode (character) VSEARCH clustering algorithm to use. One of:
#'   \itemize{
#'     \item \code{"cluster_fast"}
#'     \item \code{"cluster_size"}
#'     \item \code{"cluster_smallmem"}
#'     \item \code{"cluster_unoise"}
#'   }
#'   The selected mode is automatically converted into the corresponding flag
#'   (e.g., \code{--cluster_fast}).
#'
#' @param n_threads (integer or NULL) number of CPU threads to use, passed to
#'   VSEARCH via \code{--threads}. If \code{NULL}, threading is not explicitly
#'   specified.
#' @param minsize Minimum abundance retained by \code{cluster_unoise}. VSEARCH
#'   uses 8 by default.
#' @param unoise_alpha Abundance-skew parameter used only by
#'   \code{cluster_unoise}. VSEARCH uses 2.0 by default.
#'
#' @details
#' Mode \code{"cluster_unoise"} runs VSEARCH's UNOISE3 implementation. Input
#' labels must contain \code{size=} abundance annotations; these are generated
#' automatically for \code{sequence_table} input. UNOISE denoising does not perform
#' de novo chimera removal. Run \code{vsearch_remove_chimeras()} afterwards if
#' that step is required. The \code{identity} argument applies to the three
#' conventional clustering modes and is not passed to \code{cluster_unoise}.
#'
#' @return
#' When \code{fasta_file} is supplied, returns the exit status from
#' \code{vsearch_cmd()}. When \code{sequence_table} is supplied, returns a collapsed
#' tibble containing one row per experiment and centroid, with numeric
#' abundance columns summed across cluster members.
#'
#' @examples
#' \dontrun{
#' # Cluster the FASTA produced by chimera removal.
#' vsearch_cluster(fasta_file = "nonchimera.fasta",
#'                 output_path = "clusters", identity = 0.97, n_threads = 4)
#' }#'
#' @importFrom dplyr all_of where everything
#' @importFrom stats na.omit
#' @export
vsearch_cluster <- function(vsearch_exe_path = NULL,
                            sequence_table = NULL,
                            fasta_file = NULL,
                            output_path = NULL,
                            identity = 0.97,
                            mode = c("cluster_fast", "cluster_size", "cluster_smallmem","cluster_unoise"),
                            n_threads = NULL,
                            minsize = 8L,
                            unoise_alpha = 2.0){
  .check_number(identity, "identity", 0, 1)
  .check_number(minsize, "minsize", 1, .Machine$integer.max, TRUE)
  .check_number(unoise_alpha, "unoise_alpha", 0)
  if (!is.null(n_threads)) .check_number(n_threads, "n_threads", 1, .Machine$integer.max, TRUE)
  if (is.null(fasta_file) == is.null(sequence_table)) {
    stop("Provide either fasta_file or sequence_table.", call. = FALSE)
  }
  if (is.null(output_path)) stop("output_path must be provided.", call. = FALSE)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  mode <- match.arg(mode)

  if (!is.null(sequence_table)) {
    dere_table <- sequence_table %>%
      group_by(SEQ_ID, SEQ) %>%
      summarise(total_count = sum(total_count), .groups = "drop")
    fasta_file <- file.path(output_path, "vsearch_cluster.fasta")
    temp_seqs <- DNAStringSet(dere_table$SEQ)
    names(temp_seqs) <- paste0(dere_table$SEQ_ID, ";size=", dere_table$total_count)
    writeXStringSet(temp_seqs, fasta_file)
  }

  args <- c(paste0("--", mode), shQuote(fasta_file),
            "--centroids", shQuote(file.path(output_path, "vsearch_cluster.fas")),
            if (mode == "cluster_unoise") {
              c("--minsize", minsize, "--unoise_alpha", unoise_alpha)
            } else {
              c("--id", identity)
            },
            "--uc", shQuote(file.path(output_path, "vsearch_cluster.uc")),
            if (!is.null(n_threads)) c("--threads", n_threads))
  status <- vsearch_cmd(vsearch_exe_path = vsearch_exe_path, args = args)

  if (status != 0L) stop("VSEARCH failed with exit status ", status, ".", call. = FALSE)

  if (is.null(sequence_table)) return(status)

  temp_ucfile <- file.path(output_path, "vsearch_cluster.uc") %>%
    read_tsv(show_col_types = F,
           col_names = c("type", "cluster", "c_length","similarity","orientation","N1","N2","N3","query_ID","centoid_ID")) %>%
    mutate(query_ID = query_ID %>% str_remove(";size=.*"),
           centoid_ID = centoid_ID %>% str_remove(";size=.*"))
  H_seq <- temp_ucfile %>% filter(type == "H") %>% select(query_ID,centoid_ID) %>%
    mutate(across(all_of(everything()), ~str_remove(.x, ";size.*")))

  sequence_table %>% left_join(H_seq, by = join_by(SEQ_ID == query_ID)) %>%
    mutate(SEQ_ID = case_when(is.na(centoid_ID) ~ SEQ_ID,
                        !is.na(centoid_ID) ~ centoid_ID,
                        .default = NA),
         SEQ = case_when(
           !is.na(centoid_ID) ~ NA,
           .default = SEQ)) %>%
    group_by(SEQ_ID) %>%
    mutate(SEQ = SEQ %>% na.omit() %>% unique) %>%
    select(-centoid_ID) %>%
    # clustering
    group_by(experiment, SEQ_ID ) %>%
    summarise(SEQ = unique(SEQ), across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
}

# BLAST -------------------------------------------------
#' Run an NCBI BLAST+ command-line program
#'
#' \code{blast_cmd()} runs an NCBI BLAST+ executable with the supplied arguments
#' using \code{system2()}. The executable path and command are printed before
#' execution.
#'
#' @param base_cmd (character) name of the BLAST+ subcommand to run,
#'   such as \code{"blastn"}, \code{"blastp"}, \code{"makeblastdb"}, or
#'   \code{"tblastx"}. The executable is resolved using \code{detect_cmd()}.
#'
#' @param args (character vector or NULL) additional command-line arguments
#'   passed directly to BLAST+. Each element is treated as a separate argument.
#'   For example:
#'   \code{c("-query", "queries.fasta",
#'   "-db", "nt", "-out", "analysis_results/blast_results.tsv")}.
#'
#' @param blast_exe_path (character or NULL) optional directory or full path
#'   where BLAST+ executables are located. If \code{NULL}, the function searches
#'   for the specified \code{base_cmd} in the system PATH.
#'
#' @details
#' The function:
#' \itemize{
#'   \item detects the correct BLAST+ executable,
#'   \item prints informative messages including the resolved path,
#'   \item constructs the full BLAST+ command for reproducibility,
#'   \item executes the command using \code{system2()},
#'   \item issues a warning if BLAST exits with a non-zero status.
#' }
#' This wrapper facilitates reproducible and traceable BLAST calls in
#' automated pipelines.
#'
#' @return
#' Invisibly returns the exit status of the BLAST process (integer). A value of
#' \code{0} indicates successful execution.
#'
#' @importFrom stringr str_detect
#'
#' @examples
#' \dontrun{
#' blast_cmd(base_cmd = "blastn", args = "-version")
#' # reference.fasta contains the chosen reference sequences.
#' blast_cmd("makeblastdb", c("-in", shQuote("reference.fasta"),
#'                            "-dbtype", "nucl", "-out", shQuote("reference")))
#' }#'
#' @export
blast_cmd <- function(
    base_cmd = "blastn",
    args = NULL,
    blast_exe_path = NULL
) {
  # ---- Detect executable ----

  show_time("Searching path...")
  exe <-  detect_cmd(base_cmd,path = blast_exe_path)

  show_time(paste0("Executable file found in \"", exe, "\""))

  # ---- Run BLAST ----
  show_time("Running BLAST...")
  show_time(paste("COMMAND:", exe, paste(args, collapse = " "), "\n"))
  res <- system2(command = exe, args = args)

  if (res != 0)
    warning(base_cmd, " exited with non-zero status: ", res)
  invisible(res)
}

## BLASTN -------------------------------
#' Perform BLASTN sequence queries against a nucleotide database
#'
#' \code{blastn_cmd()} searches a nucleotide database with a query FASTA file.
#' It runs BLASTN through \code{blast_cmd()} and writes tabular hits to
#' \code{output}.
#'
#' @param fasta_file (character) path to the query FASTA file. The file should
#'   contain one or more nucleotide sequences. A complete absolute path is
#'   recommended for reproducibility.
#'
#' @param output (character) path to the output file where BLAST tabular
#'   results will be written. The function writes BLAST format 6 with additional
#'   fields (\code{qcovs}, \code{staxids}).
#' @param database (character) path to the BLAST nucleotide database
#'   (i.e., files with extensions such as \code{.nhr}, \code{.nin}, \code{.nsq}).
#'   The database must be prepared beforehand using \code{makeblastdb}.
#'
#' @param blast_exe_path (character or NULL) directory containing BLAST+
#'   executables. If \code{NULL}, \code{blast_cmd()} attempts to locate
#'   \code{blastn} via \code{detect_cmd()}.
#'
#' @param max_target_seqs (integer) maximum number of target hits to report per
#'   query sequence. Passed to BLAST as \code{-max_target_seqs}.
#'
#' @param perc_identity (numeric) minimum percent identity threshold for
#'   reporting alignments (\code{-perc_identity}).
#'
#' @param qcov_hsp_perc (numeric) minimum percentage query coverage per HSP
#'   (\code{-qcov_hsp_perc}), typically used to restrict low-coverage hits.
#'
#' @param n_threads (integer or NULL) number of CPU threads allocated to BLAST
#'   (\code{-num_threads}). If \code{NULL}, the function defaults to
#'   \code{availableCores()/4}.
#'
#' @details
#' Output format is fixed to:
#' \preformatted{
#' "6 std qcovs staxids"
#' }
#' which includes all standard BLAST columns plus query coverage (\code{qcovs})
#' and taxonomy identifiers (\code{staxids}). This makes the function suitable
#' for downstream taxonomic assignment workflows.
#'
#' @return
#' Invisibly returns the exit status of the underlying BLASTN call (integer),
#' as returned by \code{blast_cmd()}. The main result is the tabular BLAST output
#' written to the file specified by \code{output}.
#'
#' @examples
#' \dontrun{
#' # queries.fasta contains the sequences to classify; reference is a BLAST prefix.
#' blastn_cmd(fasta_file = "queries.fasta", output = "blast.tsv",
#'            database = "reference", n_threads = 4)
#' }#'
#' @export
blastn_cmd <- function(fasta_file = NULL, output = NULL, database, blast_exe_path = NULL,
                       max_target_seqs = 10, perc_identity = 90, qcov_hsp_perc = 90,
                       n_threads = NULL){
  if (is.null(fasta_file)) stop("fasta_file must be provided.", call. = FALSE)
  if (is.null(output)) stop("output or output_file must be provided.", call. = FALSE)
  show_time("Starting BLASTN")
  if(is.null(n_threads)){
    n_threads <- max(1L, floor(future::availableCores()/4))
    show_time(paste0("Using threads: ", n_threads))
  }

  .check_number(n_threads, "n_threads", 1, .Machine$integer.max, TRUE)
  .check_number(max_target_seqs, "max_target_seqs", 1, .Machine$integer.max, TRUE)
  .check_number(perc_identity, "perc_identity", 0, 100)
  .check_number(qcov_hsp_perc, "qcov_hsp_perc", 0, 100)
  .check_path(fasta_file, "fasta_file")
  .check_path(output, "output")
  .check_path(database, "database")
  status <- blast_cmd(
    base_cmd = "blastn",
    args = c(
      "-query", shQuote(fasta_file),
      "-db", shQuote(database),
      # table format
      "-outfmt", "\"6 std qcovs staxids\"",
      "-max_target_seqs", max_target_seqs,
      "-out", shQuote(output),
      "-perc_identity", perc_identity,
      "-qcov_hsp_perc", qcov_hsp_perc,
      "-num_threads",n_threads
    ),
    blast_exe_path = blast_exe_path
  )
  if (status != 0L) stop("BLASTN failed with exit status ", status, ".", call. = FALSE)
  invisible(status)
}
