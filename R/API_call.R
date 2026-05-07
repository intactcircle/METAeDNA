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
#'   \code{c("--threads", "8", "--cluster_fast", "seqs.fasta")}.
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
#'     "--derep_fulllength", "input.fasta",
#'     "--output", "output.fasta",
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
#'   input_fasta  = "seqs.fasta",
#'   output_fasta = "seqs_nonchimeras.fasta",
#'   mode         = "uchime3_denovo",
#'   output_mode  = "nonchimeras",
#'   threads      = 8
#' )
#'
#' # Reference-based chimera removal
#' vsearch_remove_chimeras(
#'   input_fasta  = "seqs.fasta",
#'   output_fasta = "seqs_filtered.fasta",
#'   mode         = "uchime_ref",
#'   ref_fasta    = "silva.fasta",
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
 # browser()
  mode <- match.arg(mode) %>% paste0("--",.)
  output_mode <- match.arg(output_mode)%>% paste0("--",.)
  if((mode ==  "uchime_ref") & is.null(ref_fasta)){stop("In uchime_ref mode, a specific reference database (ref_fasta) must be provided as input.")}
  if(is.null(ref_fasta)){
    args <- c(mode,  input_fasta, output_mode,
              output_fasta, ifelse(is.null(threads), "", paste("--threads", threads)))
  }else{
    args <- c(mode,  input_fasta, output_mode,
              output_fasta, "--db", ref_fasta ,ifelse(is.null(threads), "", paste("--threads", threads)))
  }
  vsearch_cmd( vsearch_exe_path = vsearch_exe_path,args = args)
}

## vsearch_cluster ----------------------------------
#' Cluster sequences using VSEARCH clustering algorithms
#'
#' This function is a wrapper for several VSEARCH clustering algorithms,
#' including \code{cluster_fast}, \code{cluster_size}, \code{cluster_smallmem},
#' and \code{cluster_unoise}. It constructs a VSEARCH command, executes it via
#' \code{vsearch_cmd()}, and writes cluster centroids and a \code{.uc} cluster
#' mapping file to the specified output directory.
#'
#' @param vsearch_exe_path (character or NULL) optional path to the VSEARCH
#'   executable. If \code{NULL}, the function searches for \code{"vsearch"}
#'   via \code{detect_cmd()}.
#'
#' @param input_fasta (character) path to the FASTA file containing the input
#'   sequences to be clustered.
#'
#' @param output_path (character) directory where clustering output files will
#'   be written. The function creates:
#'   \itemize{
#'     \item \code{vsearch_cluster.fas} — cluster centroid sequences,
#'     \item \code{vsearch_cluster.uc}  — VSEARCH cluster mapping file.
#'   }
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
#'
#' @return
#' Returns the exit status of the VSEARCH command executed by
#' \code{vsearch_cmd()}. The function also prints the resolved executable path
#' and the full command string.
#'
#' @examples
#' \dontrun{
#' vsearch_cluster(
#'   vsearch_exe_path = "/usr/local/bin/",
#'   input_fasta      = "seqs.fasta",
#'   output_path      = "results/",
#'   identity         = 0.97,
#'   mode             = "cluster_fast",
#'   threads          = 8
#' )
#' }
#'
#' @export
vsearch_cluster <- function(vsearch_exe_path = NULL,
                            input_fasta = NULL,
                            output_path = NULL,
                            identity = 0.97,
                            mode = c("cluster_fast", "cluster_size", "cluster_smallmem","cluster_unoise"),
                            threads = NULL){
  mode <- match.arg(mode) %>% paste0("--",.)
  temp_string <- input_fasta %>% str_remove(".*/") %>% str_remove("[.]fas.*")
  args <- c(mode,  input_fasta, "--centroids" , file.path(output_path, "vsearch_cluster.fas"), "--id", identity,
            "--uc",  file.path(output_path, "vsearch_cluster.uc"),ifelse(is.null(threads), "", paste("--threads", threads)))
  vsearch_cmd( vsearch_exe_path = vsearch_exe_path,args = args)
}

## vsearch_cluster_table--------------------------------------
#' Cluster sequences in a sequence table using VSEARCH and collapse counts accordingly
#'
#' This function performs sequence clustering directly from a sequence table
#' (as produced by \code{fasta_to_table()} or similar METAeDNA workflows).
#' Sequences are dereplicated by \code{SEQ_ID}, written to a temporary FASTA
#' file with abundance encoded in \code{size=} fields, clustered using VSEARCH,
#' and merged back into the original table according to the VSEARCH
#' \code{.uc} mapping file. The final output is a collapsed sequence table in
#' which all sequences belonging to the same centroid share a unified
#' \code{SEQ_ID} with summed read counts.
#'
#' @param seq_table (data.frame or tibble) a sequence-abundance table containing
#'   at least the columns \code{SEQ_ID}, \code{SEQ}, and \code{total_count},
#'   usually generated by earlier denoising or dereplication steps.
#'
#' @param vsearch_exe_path (character or NULL) optional path to the VSEARCH
#'   executable. If \code{NULL}, \code{vsearch_cmd()} searches for \code{"vsearch"}
#'   in the system PATH.
#'
#' @param output_path (character) directory where temporary FASTA files and the
#'   resulting VSEARCH output (\code{vsearch_cluster.fas} and
#'   \code{vsearch_cluster.uc}) will be written.
#'
#' @param identity (numeric) sequence identity threshold used for clustering.
#'   Passed to VSEARCH as \code{--id}. Typical values range from 0.90 to 0.99.
#'
#' @param mode (character) VSEARCH clustering algorithm. One of:
#'   \itemize{
#'     \item \code{"cluster_fast"}
#'     \item \code{"cluster_size"}
#'     \item \code{"cluster_smallmem"}
#'     \item \code{"cluster_unoise"}
#'   }
#'   The selected algorithm is translated into the VSEARCH flag
#'   (e.g., \code{--cluster_fast}).
#'
#' @param threads (integer or NULL) number of CPU threads to pass to VSEARCH via
#'   \code{--threads}. If \code{NULL}, threading is not explicitly set.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Dereplicates sequences by \code{SEQ_ID} and sums \code{total_count}.
#'   \item Writes a temporary FASTA file where each header includes a
#'         \code{size=} field used by VSEARCH for abundance-aware clustering.
#'   \item Executes VSEARCH clustering via \code{vsearch_cluster()}.
#'   \item Reads and parses the VSEARCH \code{.uc} mapping file to identify
#'         relationships between query sequences and their centroid clusters.
#'   \item Updates the original \code{seq_table} so that:
#'         \itemize{
#'           \item all sequences belonging to the same centroid adopt the centroid’s \code{SEQ_ID},
#'           \item sequence abundances across experiments and samples are summed,
#'           \item sequence strings (\code{SEQ}) are retained only for centroids.
#'         }
#'   \item Returns a collapsed sequence table with unified cluster membership.
#' }
#'
#' @return
#' A tibble where:
#' \itemize{
#'   \item each row corresponds to a centroid sequence cluster,
#'   \item \code{SEQ_ID} refers to the centroid identifier,
#'   \item \code{SEQ} contains the centroid sequence,
#'   \item numeric columns contain summed abundance values across all sequences
#'         assigned to the cluster,
#'   \item \code{experiment} and other metadata columns are preserved and updated.
#' }
#'
#' @examples
#' \dontrun{
#' clustered_table <- vsearch_cluster_table(
#'   seq_table       = my_seq_table,
#'   vsearch_exe_path = "/usr/local/bin",
#'   output_path      = "cluster_results/",
#'   identity         = 0.97,
#'   mode             = "cluster_fast",
#'   threads          = 8
#' )
#' }
#' @importFrom tidyselect all_of any_of
#' @export
vsearch_cluster_table <- function(seq_table,
                                  vsearch_exe_path = NULL,
                                  output_path = NULL,
                                  identity = 0.97,
                                  mode = c("cluster_fast", "cluster_size", "cluster_smallmem","cluster_unoise"),
                                  threads = NULL){
  mode <- match.arg(mode)
  # browser()

  dere_table <- seq_table %>% group_by(SEQ_ID, SEQ) %>% summarise(total_count = sum(total_count), .groups	 = "drop")
  temp_fasta_file <- file.path(output_path,"vsearch_cluster.fasta")
  temp_seqs <- dere_table$SEQ %>% DNAStringSet()
  names(temp_seqs) <- dere_table$SEQ_ID %>% paste0(";size=",dere_table$total_count)
  temp_seqs %>% writeXStringSet(temp_fasta_file)

  vsearch_cluster(vsearch_exe_path = vsearch_exe_path, input_fasta = temp_fasta_file, output_path = output_path, identity = identity, mode = mode, threads = threads)

  temp_ucfile <- file.path(output_path, "vsearch_cluster.uc") %>%
    read_tsv(show_col_types = F,
           col_names = c("type", "cluster", "c_length","similarity","orientation","N1","N2","N3","query_ID","centoid_ID")) %>%
    mutate(query_ID = query_ID %>% str_remove(";size=.*"),
           centoid_ID = centoid_ID %>% str_remove(";size=.*"))
  H_seq <- temp_ucfile %>% filter(type == "H") %>% select(query_ID,centoid_ID) %>%
    mutate(across(all_of(everything()), ~str_remove(.x, ";size.*")))

  seq_table_cluster <- seq_table %>% left_join(H_seq, by = join_by(SEQ_ID == query_ID)) %>%
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
  seq_table_cluster %>% return()
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
#'   \code{c("-query", "input.fasta", "-db", "nt", "-out", "result.txt")}.
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
#'   args = c("-query", "seqs.fasta",
#'            "-db", "nt",
#'            "-evalue", "1e-10",
#'            "-outfmt", "6",
#'            "-out", "blast_output.tsv")
#' )
#'
#' # Example: specify BLAST+ installation directory
#' blast_cmd(
#'   base_cmd = "makeblastdb",
#'   args = c("-in", "ref.fasta", "-dbtype", "nucl"),
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
#'   input   = "sequence.fasta",
#'   output  = "query_results.tsv",
#'   database = "./ITS_RefSeq_Fungi"
#' )
#'
#' # Specify BLAST installation directory
#' blastn_cmd(
#'   input   = "seqs.fasta",
#'   output  = "blast_out.tsv",
#'   database = "nt_db",
#'   blast_exe_path = "/usr/local/ncbi/blast/bin"
#' )
#' }
#'
#' @export
blastn_cmd <- function(input, output, database, blast_exe_path,max_target_seqs = 10, perc_identity = 90, qcov_hsp_perc = 90, threads = NULL){
  show_time("Starting BLASTN")
  if(is.null(threads)){
    threads <- future::availableCores()/4
    show_time(paste0("Using threads: ", threads))
  }

  message("[IMPORTANT]: Please ensure that the input path is a complete absolute path.")
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



