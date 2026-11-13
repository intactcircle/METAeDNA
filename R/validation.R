.check_number <- function(x, name, minimum = -Inf, maximum = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x > maximum || (integer && x != floor(x))) {
    stop(name, " must be one finite ", if (integer) "integer" else "number",
         " between ", minimum, " and ", maximum, ".", call. = FALSE)
  }
}

.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
}

.check_columns <- function(x, columns, name) {
  if (!is.data.frame(x) || anyDuplicated(names(x))) {
    stop(name, " must be a data frame with unique column names.", call. = FALSE)
  }
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
}

.check_counts <- function(x) {
  valid <- vapply(x, function(y) is.numeric(y) && all(is.finite(y) & y >= 0), logical(1))
  if (!all(valid)) stop("Counts must be finite non-negative numbers: ",
                        paste(names(x)[!valid], collapse = ", "), call. = FALSE)
}

.check_path <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(name, " must be one non-empty file path.", call. = FALSE)
  }
}

.check_processing <- function(arguments) {
  for (name in names(arguments)) {
    value <- arguments[[name]]
    if (name %in% c("fastq_1", "fastq_2", "fastq_file", "fasta_file", "reaction_table",
                    "out_file", "output_fasta", "output_fastq", "identified_output", "unidentified_output")) {
      .check_path(value, name)
    } else if (name %in% c("compress_output", "with_tag", "keep_internal", "keep_previous_stat", "enable_avx2", "debug")) {
      .check_flag(value, name)
    } else if (name %in% c("alpha", "beta", "threshold")) {
      .check_number(value, name)
    } else {
      .check_number(value, name, if (name %in% c("compress_level", "max_mismatches", "primer_mismatch", "max_diff", "max_indel")) 0 else 1,
                    if (name == "compress_level") 9 else .Machine$integer.max, TRUE)
    }
  }
}

# With one experiment, older control tables can omit the experiment column.
.reaction_design <- function(seq_table, control_table) {
  required <- c("group", "sample_id", "reaction_id", "category")
  .check_columns(seq_table, character(), "seq_table")
  .check_columns(control_table, required, "control_table")
  controls <- as.data.frame(control_table, stringsAsFactors = FALSE)
  if (!"experiment" %in% names(seq_table) && "experiment" %in% names(controls) &&
      length(unique(controls$experiment)) > 1L) stop("seq_table needs experiment for multiple experiments.", call. = FALSE)
  default_experiment <- if ("experiment" %in% names(controls) && nrow(controls)) as.character(controls$experiment[1]) else "experiment"
  experiment <- if ("experiment" %in% names(seq_table)) as.character(seq_table$experiment) else rep(default_experiment, nrow(seq_table))
  if (anyNA(experiment) || any(!nzchar(experiment))) stop("experiment cannot be missing or empty.", call. = FALSE)
  if (!"experiment" %in% names(controls)) {
    if (length(unique(experiment)) > 1L) stop("control_table needs experiment when seq_table contains multiple experiments.", call. = FALSE)
    controls$experiment <- if (length(experiment)) experiment[1] else "experiment"
  }
  for (column in c(required, "experiment")) {
    controls[[column]] <- trimws(as.character(controls[[column]]))
    if (anyNA(controls[[column]]) || any(!nzchar(controls[[column]]))) {
      stop("control_table$", column, " cannot be missing or empty.", call. = FALSE)
    }
  }
  controls <- unique(controls)
  if (anyDuplicated(controls[c("experiment", "reaction_id")])) {
    stop("Each experiment/reaction_id pair must identify one reaction in control_table.", call. = FALSE)
  }
  if (length(setdiff(unique(experiment), controls$experiment))) {
    stop("Some experiments are absent from control_table.", call. = FALSE)
  }
  list(table = controls, experiment = experiment)
}
