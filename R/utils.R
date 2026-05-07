# time display ----------------------------------------------
#' Print timestamped or dynamic progress messages
#'
#' This utility function provides a unified interface for printing timestamped
#' messages, updating progress in place, or displaying formatted elapsed time.
#' It is intended for use throughout the METAeDNA workflow to provide consistent
#' console output in long-running tasks.
#'
#' @param s (character) message string to print. Used in modes \code{0} and
#'   \code{2}. Ignored when \code{mode = 1}.
#'
#' @param mode (integer) controls the output format:
#'   \itemize{
#'     \item \code{0} — print a normal message with a timestamp prefix.
#'     \item \code{1} — print an elapsed time summary using the value of
#'           \code{secs}.
#'     \item \code{2} — dynamically update a single console line (typically used
#'           for streaming progress updates).
#'   }
#'
#' @param secs (numeric or NULL) elapsed time in seconds. Required only when
#'   \code{mode = 1}. Determines whether time is displayed in seconds, minutes,
#'   hours, or days.
#'
#' @details
#' \strong{Mode 0: Timestamped message}
#' Prints a message in the form
#' \preformatted{[YYYY-MM-DD HH:MM:SS] message}
#'
#' \strong{Mode 1: Elapsed time}
#' Displays the elapsed time using appropriate units:
#' seconds (< 60), minutes (< 3600), hours (< 86400), or days.
#'
#' \strong{Mode 2: Dynamic updating}
#' Updates the current console line using \code{\\r} and flushes the output.
#' Useful for real-time progress updates during large file streaming.
#'
#' @return
#' Invisibly returns \code{NULL}. The primary effect is console printing.
#'
#' @importFrom utils flush.console
#'
#' @examples
#' \dontrun{
#' show_time("Starting process...")
#'
#' # dynamic progress
#' for (i in 1:5) {
#'   show_time(paste("Progress", i, "/ 5"), mode = 2)
#'   Sys.sleep(0.5)
#' }
#' cat("\n")
#'
#' # print elapsed time
#' show_time("", mode = 1, secs = 125)
#' }
#'
show_time <- function(s, mode = 0, secs = NULL) {
  # current time
  now <- Sys.time()
  time_str <- format(now, "%Y-%m-%d %H:%M:%S")

  if (mode == 0) {
    # normal print with timestamp
    cat(sprintf("[%s] %s\n", time_str, s))

  } else if (mode == 1) {
    # elapsed time display
    if (is.null(secs)) stop("secs must be provided for mode = 1")

    cat("Elapsed time: ")

    if (secs < 60) {
      cat(sprintf("%.2f seconds", secs))
    } else if (secs < 3600) {
      cat(sprintf("%.2f minutes", secs / 60))
    } else if (secs < 86400) {
      cat(sprintf("%.2f hours", secs / 3600))
    } else {
      cat(sprintf("%.2f days", secs / 86400))
    }

    cat("\n----------------------------------------------------------------\n")

  } else if (mode == 2) {
    # dynamic one-line updating
    cat(sprintf("\r[%s] %s", time_str, s))
    flush.console()

  }
}
# find exe file path -----------------------------
#' Detect the executable path of a system command across platforms
#'
#' This utility function attempts to locate an external command-line executable
#' (e.g., \code{blastn}, \code{vsearch}) either from the system PATH or from a
#' user-specified directory. It provides cross-platform support for Windows,
#' macOS, and Linux and is used internally by wrapper functions such as
#' \code{blast_cmd()} and \code{vsearch_cmd()}.
#'
#' @param cmd (character) the name of the executable to search for (e.g.,
#'   \code{"blastn"}, \code{"vsearch"}, \code{"makeblastdb"}). Do not include
#'   an extension such as \code{.exe}.
#'
#' @param path (character or NULL) optional directory in which to search for the
#'   executable. If \code{NULL}, the function first queries the system PATH via
#'   \code{Sys.which()} and then searches several common installation locations
#'   depending on the operating system.
#'
#' @details
#' The detection logic proceeds in the following order:
#' \enumerate{
#'   \item Use \code{Sys.which(cmd)} to check for the executable in the
#'         system PATH.
#'   \item If not found, search platform-specific default directories:
#'         \itemize{
#'           \item Windows: typical BLAST+/VSEARCH installation folders under
#'                 \code{"C:/Program Files"} and \code{"C:/blast"}.
#'           \item macOS: \code{/usr/local/ncbi/blast/bin},
#'                 \code{/opt/homebrew/bin},
#'                 \code{/usr/local/bin}.
#'           \item Linux: \code{/usr/local/bin}, \code{/usr/bin},
#'                 \code{/opt/ncbi/blast+/bin}.
#'         }
#'   \item If \code{path} is supplied, construct the executable path using the
#'         directory provided by the user.
#' }
#'
#' On Windows, the function automatically appends \code{.exe} when constructing
#' candidate executable paths.
#'
#' If the executable cannot be located, the function throws an error with a
#' user-readable message instructing how to provide the correct path.
#'
#' @return
#' Returns a character string giving the full path to the located executable.
#' An error is thrown if the command cannot be found.
#'
#' @importFrom stringr str_detect
#'
#' @examples
#' \dontrun{
#' # Detect vsearch from system PATH
#' detect_cmd("vsearch")
#'
#' # Detect blastn from a user-specified installation directory
#' detect_cmd("blastn", path = "/usr/local/ncbi/blast/bin")
#' }
#'
detect_cmd <- function(cmd, path = NULL) {
  sys <- Sys.info()[["sysname"]]
  if(is.null(path)){
    found <- suppressWarnings(Sys.which(cmd))
    if (nzchar(found)) return(found)
    paths <- switch(sys,
                    "Windows" = c("C:/Program Files/NCBI/blast+/bin", "C:/blast+/bin", "C:/Program Files/NCBI/blast+/bin"),
                    "Darwin"  = c("/usr/local/ncbi/blast/bin", "/opt/homebrew/bin", "/usr/local/bin"),
                    "Linux"   = c("/usr/local/bin", "/usr/bin", "/opt/ncbi/blast+/bin"),
                    NULL
    )
    for (p in paths) {
      if(str_detect(sys, "Windows|windows") ){
        full <-paste0(p,"/", cmd, ".exe")
      }else{
        full <- file.path(p, cmd)
      }
      if (file.exists(full)) return(full)
    }
    stop(paste0("Cannot find [", cmd,"] command.\nPlease ensure that program are installed and correctly configured within the system, or manually specify the location of the executable file." ), call. = FALSE)

  }else{
    if(str_detect(sys, "Windows|windows") ){
      full <-paste0(path,"/", cmd, ".exe")
    }else{
      full <- file.path(path, cmd)
    }
    if(file.exists(full)){return(full)}else{ stop(paste0("Cannot find [", cmd,"] command.\nPlease ensure specify a correct location of the executable file." ), call. = FALSE)
}
  }
}
