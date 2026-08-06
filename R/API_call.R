# Vsearch -------------------------------------
## vsearch command --------------------------------
#' Execute a VSEARCH command with automatic executable detection
#'
#' This function locates the VSEARCH executable (optionally using a user-supplied
#' path), constructs a command, prints the resolved executable location and full
#' command string, and runs it using \code{system2()}. It is intended as a simple
#' wrapper for running VSEARCH subcommands (e.g., \code{--cluster_fast},
#' \code{--uchime3_denovo}, \code{--search_global}, etc.) within the METAeDNA
#' workflow.
#'
#' @param vsearch_exe_path (character or NULL) optional path to the VSEARCH
#'   executable. If \code{NULL}, the function searches for \code{"vsearch"} in
#'   the system PATH using \code{detect_cmd()}. If a directory is provided, the
#'   function attempts to locate the binary inside that directory.
#'
#' @param args (character vector) command-line arguments passed directly to
#'   \code{system2()}. Each element is treated as one argument; for example:
#'   \code{c("--threads", "8", "--cluster_fast",
#'   "example_data/sequences.fasta")}.
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
#' # Example: run VSEARCH dereplication
#' vsearch_cmd(
#'   args = c(
#'     "--derep_fulllength", "example_data/sequences.fasta",
#'     "--output", "analysis_results/unique_sequences.fasta",
#'     "--threads", "4"
#'   )
#' )
#'
#' # Example: explicitly specify executable path
#' vsearch_cmd(
#'   vsearch_exe_path = "/usr/local/bin",
#'   args = c("--version")
#' )
#' }
#'
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
#' This function provides a wrapper for several VSEARCH chimera detection
#' algorithms, including \code{uchime_denovo}, \code{uchime2_denovo},
#' \code{uchime3_denovo}, and \code{uchime_ref}. It automatically constructs and
#' executes the corresponding VSEARCH command, optionally using a user-specified
#' VSEARCH executable path. Both reference-free and reference-based chimera
#' detection modes are supported.
#'
#' @param vsearch_exe_path (character or NULL) optional path to the VSEARCH
#'   executable. If \code{NULL}, the function searches for \code{"vsearch"} on
#'   the system PATH using \code{detect_cmd()}.
#'
#' @param input_fasta (character) path to the input FASTA file containing the
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
#' @param threads (integer or NULL) number of CPU threads to pass to VSEARCH via
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
#' # De novo chimera removal (uchime3_denovo)
#' vsearch_remove_chimeras(
#'   input_fasta  = "example_data/sequences.fasta",
#'   output_fasta = "analysis_results/nonchimeric_sequences.fasta",
#'   mode         = "uchime3_denovo",
#'   output_mode  = "nonchimeras",
#'   threads      = 8
#' )
#'
#' # Reference-based chimera removal
#' vsearch_remove_chimeras(
#'   input_fasta  = "example_data/sequences.fasta",
#'   output_fasta = "analysis_results/filtered_sequences.fasta",
#'   mode         = "uchime_ref",
#'   ref_fasta    = "reference_databases/silva.fasta",
#'   output_mode  = "nonchimeras"
#' )
#' }
#'
#' @export
vsearch_remove_chimeras <- function(
    vsearch_exe_path = NULL,
    input_fasta,
    output_fasta,
    mode = c("uchime_denovo", "uchime2_denovo", "uchime3_denovo", "uchime_ref"),
    output_mode = c("nonchimeras", "chimeras" ,   "uchimealns", "uchimeout"),
    ref_fasta = NULL,
    minuniquesize = NULL,
    threads = NULL
){
  mode <- match.arg(mode)
  output_mode <- match.arg(output_mode)
  if((mode ==  "uchime_ref") & is.null(ref_fasta)){stop("In uchime_ref mode, a specific reference database (ref_fasta) must be provided as input.")}
  mode <- paste0("--", mode)
  output_mode <- paste0("--", output_mode)
  if(is.null(ref_fasta)){
    args <- c(mode, input_fasta, output_mode, output_fasta,
              if (!is.null(threads)) c("--threads", threads))
  }else{
    args <- c(mode, input_fasta, output_mode, output_fasta, "--db", ref_fasta,
              if (!is.null(threads)) c("--threads", threads))
  }
  vsearch_cmd( vsearch_exe_path = vsearch_exe_path,args = args)
}

## vsearch_cluster ----------------------------------
#' Cluster a FASTA file or sequence table using VSEARCH
#'
#' The function accepts either a FASTA file through \code{input_fasta} or a
#' METAeDNA sequence-abundance table through \code{seq_table}. Both input modes
#' write cluster centroids and a \code{.uc} mapping file. When \code{seq_table}
#' is supplied, cluster membership is additionally merged back into the table
#' and numeric abundance columns are summed for each experiment and centroid.
#'
#' @param vsearch_exe_path (character or NULL) optional path to the VSEARCH
#'   executable. If \code{NULL}, the function searches for \code{"vsearch"}
#'   via \code{detect_cmd()}.
#'
#' @param input_fasta (character) path to the FASTA file containing the input
#'   sequences. Supply either \code{input_fasta} or \code{seq_table}.
#' @param seq_table (data.frame or tibble) optional sequence-abundance table
#'   containing \code{SEQ_ID}, \code{SEQ}, and \code{total_count}. Supply either
#'   \code{seq_table} or \code{input_fasta}.
#'
#' @param output_path (character) directory where clustering output files will
#'   be written. The function creates:
#'   \itemize{
#'     \item \code{vsearch_cluster.fas} — cluster centroid sequences,
#'     \item \code{vsearch_cluster.uc}  — VSEARCH cluster mapping file.
#'   }
#' @param output_dir Deprecated alias of \code{output_path}.
#'
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
#' @param threads (integer or NULL) number of CPU threads to use, passed to
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
#' automatically for \code{seq_table} input. UNOISE denoising does not perform
#' de novo chimera removal. Run \code{vsearch_remove_chimeras()} afterwards if
#' that step is required. The \code{identity} argument applies to the three
#' conventional clustering modes and is not passed to \code{cluster_unoise}.
#'
#' @return
#' When \code{input_fasta} is supplied, returns the exit status from
#' \code{vsearch_cmd()}. When \code{seq_table} is supplied, returns a collapsed
#' tibble containing one row per experiment and centroid, with numeric
#' abundance columns summed across cluster members.
#'
#' @examples
#' \dontrun{
#' vsearch_cluster(
#'   vsearch_exe_path = "/usr/local/bin/",
#'   input_fasta      = "example_data/sequences.fasta",
#'   output_path      = "analysis_results/sequence_clusters",
#'   identity         = 0.97,
#'   mode             = "cluster_fast",
#'   threads          = 8
#' )
#'
#' clustered_table <- vsearch_cluster(
#'   seq_table  = my_seq_table,
#'   output_path = "analysis_results/table_clusters",
#'   identity   = 0.97,
#'   threads    = 8
#' )
#' }
#'
#' @importFrom tidyselect all_of where
#' @export
vsearch_cluster <- function(vsearch_exe_path = NULL,
                            input_fasta = NULL,
                            output_path = NULL,
                            identity = 0.97,
                            mode = c("cluster_fast", "cluster_size", "cluster_smallmem","cluster_unoise"),
                            threads = NULL,
                            output_dir = output_path,
                            seq_table = NULL,
                            minsize = 8L,
                            unoise_alpha = 2.0){
  if (is.null(output_path)) output_path <- output_dir
  if (is.null(input_fasta) == is.null(seq_table)) {
    stop("Provide either input_fasta or seq_table.", call. = FALSE)
  }
  if (is.null(output_path)) stop("output_path or output_dir must be provided.", call. = FALSE)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  mode <- match.arg(mode)

  if (!is.null(seq_table)) {
    dere_table <- seq_table %>%
      group_by(SEQ_ID, SEQ) %>%
      summarise(total_count = sum(total_count), .groups = "drop")
    input_fasta <- file.path(output_path, "vsearch_cluster.fasta")
    temp_seqs <- DNAStringSet(dere_table$SEQ)
    names(temp_seqs) <- paste0(dere_table$SEQ_ID, ";size=", dere_table$total_count)
    writeXStringSet(temp_seqs, input_fasta)
  }

  args <- c(paste0("--", mode), input_fasta,
            "--centroids", file.path(output_path, "vsearch_cluster.fas"),
            if (mode == "cluster_unoise") {
              c("--minsize", minsize, "--unoise_alpha", unoise_alpha)
            } else {
              c("--id", identity)
            },
            "--uc", file.path(output_path, "vsearch_cluster.uc"),
            if (!is.null(threads)) c("--threads", threads))
  status <- vsearch_cmd(vsearch_exe_path = vsearch_exe_path, args = args)

  if (is.null(seq_table)) return(status)

  temp_ucfile <- file.path(output_path, "vsearch_cluster.uc") %>%
    read_tsv(show_col_types = F,
           col_names = c("type", "cluster", "c_length","similarity","orientation","N1","N2","N3","query_ID","centoid_ID")) %>%
    mutate(query_ID = query_ID %>% str_remove(";size=.*"),
           centoid_ID = centoid_ID %>% str_remove(";size=.*"))
  H_seq <- temp_ucfile %>% filter(type == "H") %>% select(query_ID,centoid_ID) %>%
    mutate(across(all_of(everything()), ~str_remove(.x, ";size.*")))

  seq_table %>% left_join(H_seq, by = join_by(SEQ_ID == query_ID)) %>%
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
#' This function wraps the execution of NCBI BLAST+ command-line tools
#' (e.g. \code{blastn}, \code{blastp}, \code{blastx}, \code{tblastn},
#' \code{makeblastdb}). It automatically detects the BLAST executable,
#' constructs a command, prints the resolved executable location and full
#' command string, and executes the command using \code{system2()}.
#'
#' @param base_cmd (character) name of the BLAST+ subcommand to run,
#'   such as \code{"blastn"}, \code{"blastp"}, \code{"makeblastdb"}, or
#'   \code{"tblastx"}. The executable is resolved using \code{detect_cmd()}.
#'
#' @param args (character vector or NULL) additional command-line arguments
#'   passed directly to BLAST+. Each element is treated as a separate argument.
#'   For example:
#'   \code{c("-query", "example_data/query_sequences.fasta",
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
#' # Example: run BLASTN alignment
#' blast_cmd(
#'   base_cmd = "blastn",
#'   args = c("-query", "example_data/query_sequences.fasta",
#'            "-db", "reference_databases/blast/nt",
#'            "-evalue", "1e-10",
#'            "-outfmt", "6",
#'            "-out", "analysis_results/blast_results.tsv")
#' )
#'
#' # Example: specify BLAST+ installation directory
#' blast_cmd(
#'   base_cmd = "makeblastdb",
#'   args = c("-in", "example_data/reference_sequences.fasta",
#'            "-dbtype", "nucl"),
#'   blast_exe_path = "/usr/local/ncbi/blast/bin"
#' )
#' }
#'
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
#' This function provides a high-level interface for running BLASTN searches
#' using NCBI BLAST+. It constructs a standardized BLASTN command, ensures
#' appropriate output formatting, prints diagnostic messages, and executes the
#' command through \code{blast_cmd()}. The resulting BLAST output is written
#' directly to a user-specified file.
#'
#' @param input (character) path to the query FASTA file. The file should
#'   contain one or more nucleotide sequences. A complete absolute path is
#'   recommended for reproducibility.
#'
#' @param output (character) path to the output file where BLAST tabular
#'   results will be written. The function writes BLAST format 6 with additional
#'   fields (\code{qcovs}, \code{staxids}).
#' @param input_fasta Deprecated alias of \code{input}.
#' @param output_file Deprecated alias of \code{output}.
#'
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
#' @param threads (integer or NULL) number of CPU threads allocated to BLAST
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
#' blastn_cmd(
#'   input   = "example_data/query_sequences.fasta",
#'   output  = "analysis_results/blast_results.tsv",
#'   database = "reference_databases/ITS_RefSeq_Fungi"
#' )
#'
#' # Specify BLAST installation directory
#' blastn_cmd(
#'   input   = "example_data/query_sequences.fasta",
#'   output  = "analysis_results/blast_results.tsv",
#'   database = "reference_databases/blast/nt",
#'   blast_exe_path = "/usr/local/ncbi/blast/bin"
#' )
#' }
#'
#' @export
blastn_cmd <- function(input = NULL, output = NULL, database, blast_exe_path = NULL,
                       max_target_seqs = 10, perc_identity = 90, qcov_hsp_perc = 90,
                       threads = NULL, input_fasta = input, output_file = output){
  if (is.null(input)) input <- input_fasta
  if (is.null(output)) output <- output_file
  if (is.null(input)) stop("input or input_fasta must be provided.", call. = FALSE)
  if (is.null(output)) stop("output or output_file must be provided.", call. = FALSE)
  show_time("Starting BLASTN")
  if(is.null(threads)){
    threads <- future::availableCores()/4
    show_time(paste0("Using threads: ", threads))
  }

  blast_cmd(
    base_cmd = "blastn",
    args = c(
      "-query", input,
      "-db", database,
      # table format
      "-outfmt", "\"6 std qcovs staxids\"",
      "-max_target_seqs", max_target_seqs,
      "-out", output,
      "-perc_identity", perc_identity,
      "-qcov_hsp_perc", qcov_hsp_perc,
      "-num_threads",threads
    ),
    blast_exe_path = blast_exe_path
  )
}
