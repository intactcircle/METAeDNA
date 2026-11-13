#' Correct read counts using blank reactions
#'
#' Subtract blank counts or exclude detections using the blank reactions
#' assigned to each sample in the control table.
#'
#' @param sequence_table Sequence table containing `SEQ` and `count_TAG` columns.
#'   Include `experiment` when more than one experiment is present. `ID` is
#'   optional. Rows, including repeated sequence identifiers, are preserved.
#' @param control_table Reaction definitions with `group`, `sample_id`, `reaction_id`,
#'   `category`, and one column per blank type listed in `control_names`.
#'   Sample reactions use `category = "sample"`. Blank columns contain the
#'   corresponding blank sample name. Include `experiment` for multiple
#'   experiments; each experiment/reaction_id pair must identify one reaction.
#' @param method `"maximum"` subtracts the largest blank-replicate count.
#'   `"mean"` subtracts the rounded mean. `"threshold"` excludes detections
#'   when any blank replicate reaches `threshold`; `"both_threshold"`
#'   excludes them when every blank replicate reaches it.
#' @param control_names Blank categories to use, such as `c("FB", "EB", "PB")`.
#' @param threshold Non-negative read-count cutoff for the threshold methods.
#' @param remove_status Remove columns beginning with `status_`.
#' @param missing_columns `"zero"` fills absent count columns referenced by
#'   the control table with zero. `"error"` reports missing columns instead.
#'
#' @details Corrections are computed within each experiment and control group.
#'   With several blank types, each correction is applied in turn using the
#'   original blank counts. Subtracted counts are bounded below by zero.
#'   A blank name without a matching control-table row is an error.
#'   Control counts are cleared and control-only columns are removed.
#' @return A sequence table with corrected counts and, if present, updated
#'   `total_count`. This method does not estimate contamination probabilities.
#' @examples
#' sequence_table <- data.frame(SEQ = c("ACGT", "TGCA"), count_A = c(20, 5),
#'                         count_B = c(3, 8), total_count = c(23, 13))
#' controls <- data.frame(group = "G1", sample_id = c("S1", "blank"),
#'                        reaction_id = c("A", "B"), category = c("sample", "FB"),
#'                        FB = c("blank", NA))
#' decontamination(sequence_table, controls, control_names = "FB")
#' @export
decontamination <- function(sequence_table, control_table,
                            method = c("maximum", "threshold", "mean", "both_threshold"),
                            control_names = NULL, threshold = 10,
                            remove_status = TRUE,
                            missing_columns = c("zero", "error")) {
  method <- match.arg(method)
  missing_columns <- match.arg(missing_columns)
  .check_number(threshold, "threshold", 0)
  .check_flag(remove_status, "remove_status")
  .check_columns(sequence_table, "SEQ", "sequence_table")
  .check_columns(control_table, c("group", "sample_id", "reaction_id", "category"), "control_table")
  if (!is.character(control_names) || !length(control_names) ||
      anyNA(control_names) || any(!nzchar(control_names)) || anyDuplicated(control_names)) {
    stop("control_names must contain distinct blank category names.", call. = FALSE)
  }
  .check_columns(control_table, control_names, "control_table")
  design <- .reaction_design(sequence_table, control_table)
  controls <- design$table
  original_counts <- sequence_table
  count_columns <- paste0("count_", unique(controls$reaction_id[controls$category %in% c("sample", control_names)]))
  missing <- setdiff(count_columns, names(sequence_table))
  if (length(missing)) {
    if (missing_columns == "error") {
      stop("Missing count column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    sequence_table[missing] <- 0
    original_counts[missing] <- 0
  }
  .check_counts(sequence_table[, startsWith(names(sequence_table), "count_"), drop = FALSE])

  for (experiment in unique(controls$experiment)) {
    rows <- which(design$experiment == experiment)
    experiment_controls <- controls[controls$experiment == experiment, , drop = FALSE]
    samples <- experiment_controls[experiment_controls$category == "sample", , drop = FALSE]
    for (blank_type in control_names) {
      blank_names <- trimws(as.character(samples[[blank_type]]))
      blank_groups <- unique(data.frame(group = samples$group, blank = blank_names))
      blank_groups <- blank_groups[!is.na(blank_groups$blank) & nzchar(blank_groups$blank), , drop = FALSE]
      for (i in seq_len(nrow(blank_groups))) {
        blank_name <- blank_groups$blank[i]
        blanks <- experiment_controls[
          experiment_controls$group == blank_groups$group[i] &
          experiment_controls$sample_id == blank_name &
          experiment_controls$category == blank_type, , drop = FALSE
        ]
        if (!nrow(blanks)) {
          stop("Blank '", blank_name, "' is absent from control_table for experiment '",
               experiment, "', group '", blank_groups$group[i], "'.", call. = FALSE)
        }
        if (!length(rows)) next
        # Reactions sharing a blank use the same correction vector.
        blank_counts <- original_counts[rows, paste0("count_", blanks$reaction_id), drop = FALSE]
        correction <- switch(method,
          maximum = Reduce(pmax, blank_counts),
          mean = round(rowMeans(as.matrix(blank_counts))),
          threshold = Reduce(pmax, blank_counts) >= threshold,
          both_threshold = Reduce(pmin, blank_counts) >= threshold
        )
        reaction_ids <- samples$reaction_id[
          samples$group == blank_groups$group[i] & !is.na(blank_names) & blank_names == blank_name
        ]
        for (column in paste0("count_", reaction_ids)) {
          if (method %in% c("maximum", "mean")) {
            sequence_table[[column]][rows] <- pmax(sequence_table[[column]][rows] - correction, 0)
          } else {
            sequence_table[[column]][rows[correction]] <- 0
          }
        }
      }
    }
  }
  blank_columns <- paste0("count_", unique(controls$reaction_id[controls$category %in% control_names]))
  # A tag may be reused by a sample in another experiment.
  sample_columns <- paste0("count_", unique(controls$reaction_id[controls$category == "sample"]))
  for (experiment in unique(controls$experiment)) {
    blank_tags <- controls$reaction_id[controls$experiment == experiment & controls$category %in% control_names]
    rows <- which(design$experiment == experiment)
    for (column in intersect(paste0("count_", blank_tags), names(sequence_table))) sequence_table[[column]][rows] <- 0
  }
  sequence_table <- sequence_table[, setdiff(names(sequence_table), setdiff(blank_columns, sample_columns)), drop = FALSE]
  if (remove_status) sequence_table <- sequence_table[, !startsWith(names(sequence_table), "status_"), drop = FALSE]
  if ("total_count" %in% names(sequence_table)) {
    sequence_table$total_count <- rowSums(sequence_table[, startsWith(names(sequence_table), "count_"), drop = FALSE])
  }
  sequence_table
}


#
#' Merge PCR replicates within each experiment and sample
#'
#' Normalize each reaction to relative read abundance, then combine reactions
#' belonging to the same experiment, group and sample.
#'
#' @param sequence_table Sequence table with `count_TAG` columns. Include
#'   `experiment` when the table contains several experiments.
#' @param control_table Reaction definitions containing `group`, `sample_id`,
#'   `reaction_id` and `category`. Only `category = "sample"` is merged.
#'   Multiple experiments require an `experiment` column in both tables.
#'   A reaction ID may be reused in different experiments. Reaction IDs name
#'   count columns; they are not the oligonucleotide tag sequences.
#' @param method `"mean_rra"` averages reaction relative abundances;
#'   `"max_rra"` takes their maximum; `"both_occur"` averages them and
#'   retains only sequences found in at least `repeat_times` reactions.
#' @param repeat_times Minimum number of positive reactions for `"both_occur"`.
#' @param replicate_number Positive integer denominator for averaging.
#'   `NULL` uses the number of reaction columns for each sample. A value
#'   of 3 divides by 3 even if the control table lists only one reaction.
#' @param present_absent Convert merged abundances to 0 or 1.
#' @param remove_all_zero Remove rows absent from all merged samples.
#'
#' @details Each reaction is normalized using reads from its own experiment.
#'   Empty reactions contribute zero. Declared reaction columns must exist.
#'   Retained values are not renormalized after occurrence filtering; column
#'   sums may differ from one, especially with a fixed denominator or max RRA.
#' @return The sequence table with sample reaction columns replaced by
#'   `group_sample` columns for one experiment, or
#'   `experiment_group_sample` for several experiments. Other-experiment rows
#'   receive zero in each sample column. Status columns are removed.
#' @examples
#' sequence_table <- data.frame(experiment = rep(c("run1", "run2"), each = 2),
#'                         SEQ = rep(c("ACGT", "TGCA"), 2),
#'                         count_A = c(9, 1, 2, 8), count_B = c(3, 7, 4, 6))
#' controls <- expand.grid(experiment = c("run1", "run2"), reaction_id = c("A", "B"),
#'                         stringsAsFactors = FALSE)
#' controls$group <- "G1"
#' controls$sample_id <- "S1"
#' controls$category <- "sample"
#' merge_replicate(sequence_table, controls)
#' merge_replicate(sequence_table, controls, replicate_number = 3)
#' @export
merge_replicate <- function(sequence_table, control_table,
                            method = c("mean_rra", "max_rra", "both_occur"),
                            repeat_times = 1, replicate_number = NULL,
                            present_absent = FALSE, remove_all_zero = TRUE) {
  method <- match.arg(method)
  .check_number(repeat_times, "repeat_times", 1, integer = TRUE)
  if (!is.null(replicate_number)) .check_number(replicate_number, "replicate_number", 1, integer = TRUE)
  .check_flag(present_absent, "present_absent")
  .check_flag(remove_all_zero, "remove_all_zero")
  design <- .reaction_design(sequence_table, control_table)
  samples <- design$table[design$table$category == "sample", , drop = FALSE]
  if (!nrow(samples)) stop("control_table has no sample reactions.", call. = FALSE)
  columns <- unique(paste0("count_", samples$reaction_id))
  .check_columns(sequence_table, columns, "sequence_table")
  .check_counts(sequence_table[columns])
  pairs <- unique(samples[c("experiment", "group", "sample_id")])
  column_names <- paste(pairs$group, pairs$sample_id, sep = "_")
  if (length(unique(pairs$experiment)) > 1L) {
    column_names <- paste(pairs$experiment, column_names, sep = "_")
  }
  if (anyDuplicated(column_names) || any(column_names %in% names(sequence_table))) {
    stop("Merged sample names are duplicated or already exist in sequence_table.", call. = FALSE)
  }
  for (i in seq_len(nrow(pairs))) {
    reaction <- samples[samples$experiment == pairs$experiment[i] &
                        samples$group == pairs$group[i] &
                        samples$sample_id == pairs$sample_id[i], , drop = FALSE]
    rows <- which(design$experiment == pairs$experiment[i])
    counts <- as.matrix(sequence_table[rows, paste0("count_", reaction$reaction_id), drop = FALSE])
    totals <- colSums(counts)
    rra <- sweep(counts, 2, ifelse(totals == 0, 1, totals), "/")
    denominator <- if (is.null(replicate_number)) ncol(counts) else replicate_number
    values <- if (method == "max_rra" && length(rows)) apply(rra, 1, max) else rowSums(rra) / denominator
    if (method == "both_occur") values[rowSums(counts > 0) < repeat_times] <- 0
    if (present_absent) values <- as.numeric(values > 0)
    sequence_table[[column_names[i]]] <- 0
    sequence_table[[column_names[i]]][rows] <- values
  }
  sequence_table <- sequence_table[, !names(sequence_table) %in% columns & !startsWith(names(sequence_table), "status_"), drop = FALSE]
  if (remove_all_zero) sequence_table <- sequence_table[rowSums(sequence_table[column_names]) > 0, , drop = FALSE]
  sequence_table
}

#' Construct a taxon abundance matrix from merged reactions
#'
#' Match classification by sequence and sum abundances assigned to the same
#' scientific name.
#'
#' @param merged_table A sequence table returned by `merge_replicate()`,
#'   containing `SEQ` and the merged sample-abundance columns.
#' @param control_table The control table used for merging replicates. Rows
#'   with `category = "sample"` identify the abundance columns through `group`
#'   and `sample_id`, with an experiment prefix when several experiments are
#'   present. Other numeric columns, such as `total_count`, are not included
#'   in the matrix.
#' @param taxonomy_table Classification table with `SEQ` and
#'   `scientific_name`. It supplies the taxon name matched to each `SEQ` in
#'   `merged_table`; abundances sharing a name are summed. Each sequence in
#'   `merged_table` must have one non-missing classification.
#' @param taxa_in_rows If `FALSE`, samples are rows and taxa are columns.
#'   Use `TRUE` to transpose this arrangement.
#' @details
#' The three tables supply sequence abundances, sample definitions, and
#' classifications, respectively. `control_table` selects the sample columns;
#' `taxonomy_table` supplies their taxon grouping through `SEQ`. Both tables
#' are required, even if `merged_table` already contains classification columns.
#' The function does not repeat blank correction, replicate merging, or
#' relative-abundance normalization. Matrix values retain the abundance scale
#' of `merged_table`.
#' @return A numeric matrix. Sample names follow `merge_replicate()`,
#'   including experiment prefixes when several experiments are present.
#' @examples
#' merged <- data.frame(SEQ = c("ACGT", "TGCA"), G1_S1 = c(0.3, 0.7))
#' controls <- data.frame(group = "G1", sample_id = "S1", reaction_id = "A",
#'                        category = "sample")
#' taxonomy <- data.frame(SEQ = merged$SEQ,
#'                        scientific_name = c("Species alpha", "Species beta"))
#' merge_to_matrix(merged, controls, taxonomy)
#' @export
merge_to_matrix <- function(merged_table, control_table, taxonomy_table,
                            taxa_in_rows = FALSE) {
  required <- c("group", "sample_id", "category")
  missing <- setdiff(required, names(control_table))
  if (length(missing)) {
    stop("control_table is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  .check_flag(taxa_in_rows, "taxa_in_rows")
  design <- .reaction_design(merged_table, control_table)
  sample_table <- design$table[design$table$category == "sample", , drop = FALSE]
  sample_columns <- paste(sample_table$group, sample_table$sample_id, sep = "_")
  if (length(unique(sample_table$experiment)) > 1L) {
    sample_columns <- paste(sample_table$experiment, sample_columns, sep = "_")
  }
  sample_columns <- unique(sample_columns)
  if (!length(sample_columns)) {
    stop("control_table does not contain any rows with category = 'sample'.",
         call. = FALSE)
  }
  missing_samples <- setdiff(sample_columns, names(merged_table))
  if (length(missing_samples)) {
    stop("merged_table is missing sample column(s): ",
         paste(missing_samples, collapse = ", "), call. = FALSE)
  }

  if (!"SEQ" %in% names(merged_table)) {
    stop("merged_table is missing: SEQ", call. = FALSE)
  }
  taxonomy_columns <- c("SEQ", "scientific_name")
  if (!all(taxonomy_columns %in% names(taxonomy_table))) {
    stop("taxonomy_table must contain SEQ and scientific_name.", call. = FALSE)
  }

  valid_samples <- vapply(merged_table[sample_columns], function(x) {
    is.numeric(x) && !anyNA(x) && all(x >= 0)
  }, logical(1))
  if (!all(valid_samples)) {
    stop("Sample columns must contain non-negative numbers without missing values: ",
         paste(sample_columns[!valid_samples], collapse = ", "), call. = FALSE)
  }

  taxonomy_table <- unique(taxonomy_table[c("SEQ", "scientific_name")])
  if (anyDuplicated(taxonomy_table$SEQ)) stop("Each SEQ must have one scientific_name.", call. = FALSE)
  scientific_name <- taxonomy_table$scientific_name[
    match(merged_table$SEQ, taxonomy_table$SEQ)
  ]
  if (anyNA(scientific_name) || any(!nzchar(scientific_name))) {
    stop("Some sequences do not have a scientific_name.",
         call. = FALSE)
  }

  abundance <- as.matrix(merged_table[, sample_columns, drop = FALSE])
  abundance <- rowsum(abundance, group = scientific_name, reorder = FALSE)
  if (taxa_in_rows) abundance else t(abundance)
}


#' Calculate read-based rarefaction curves for each reaction
#'
#' Calculate the expected number of taxa recovered at increasing sequencing
#' depths. The function can use either a taxa-by-reaction count matrix or a
#' decontaminated sequence table together with its classification result.
#'
#' @param sequence_table A numeric taxa-by-reaction matrix, or a decontaminated
#'   sequence table containing `SEQ` and `count_` columns.
#' @param taxonomy_table Optional classification result containing `SEQ` and
#'   `scientific_name`. The two input tables are matched only through `SEQ`.
#' @param n_threads Number of parallel processes. Default is 1.
#' @param bin_size Sequencing-depth interval between adjacent curve points.
#'   When `NULL`, each reaction uses approximately 100 log-spaced points so
#'   that the rapidly changing low-depth part of the curve is represented.
#' @param endpoint Optional extrapolation endpoint. By default, extrapolation
#'   for every reaction ends at `maximum + maximum * 0.1`, where `maximum` is
#'   the largest observed sequencing depth in the input table.
#' @param x_limit,y_limit Optional axis limits.
#' @param show_legend Legend position. Use `NULL` to hide the legend.
#' @param show_observed Show the observed richness at the full read depth.
#' @param show_extrapolation Show richness extrapolation beyond the observed
#'   sequencing depth as a dashed line.
#' @param x_title,y_title Axis titles.
#' @param facet_experiment Split the plot by `experiment` when the sequence
#'   table contains that column. Reused reaction IDs are calculated separately
#'   within each experiment; the default extrapolation endpoint is set within
#'   each panel. `FALSE` pools experiments. Ignored for matrix input.
#'
#' @return A list containing `rarefaction_data` and the ggplot object `fig`.
#'
#' @importFrom dplyr bind_rows
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom ggplot2 aes coord_cartesian element_rect geom_line geom_point
#'   ggplot labs scale_colour_manual scale_linetype_manual theme theme_classic
#'   facet_wrap element_blank
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' set.seed(123)
#' reaction_reads <- round(seq(250, 2000, length.out = 20))
#' species_number <- round(seq(6, 30, length.out = 20))
#' counts <- vapply(seq_len(20), function(i) {
#'   species <- species_number[i]
#'   species_reads <- as.vector(rmultinom(
#'     1, reaction_reads[i], prob = runif(species, 0.5, 1.5)
#'   ))
#'   sequence_reads <- unlist(lapply(species_reads, function(x) {
#'     first <- rbinom(1, x, 0.5)
#'     c(first, x - first)
#'   }))
#'   c(sequence_reads, rep(0, 60 - length(sequence_reads)))
#' }, numeric(60))
#' colnames(counts) <- sprintf("count_F%03d", seq_len(20))
#' sequence_table <- data.frame(
#'   SEQ = paste0("sequence_", 1:60), counts, check.names = FALSE
#' )
#' taxonomy <- data.frame(
#'   SEQ = sequence_table$SEQ,
#'   scientific_name = rep(paste("Species", 1:30), each = 2)
#' )
#' reaction_info <- data.frame(
#'   tag = sub("^count_", "", colnames(counts)),
#'   reads = reaction_reads,
#'   expected_species = species_number
#' )
#' result <- rarefaction_extrapolation(
#'   sequence_table, taxonomy, n_threads = 2, bin_size = 20,
#'   show_extrapolation = TRUE, show_legend = "right"
#' )
#' result$fig
#'
#' @export
rarefaction_extrapolation <- function(sequence_table, taxonomy_table = NULL,
                                      n_threads = 1, bin_size = NULL,
                                      endpoint = NULL, x_limit = NULL,
                                      y_limit = NULL, show_legend = NULL,
                                      show_observed = FALSE,
                                      show_extrapolation = FALSE,
                                      x_title = "Sequencing depth",
                                      y_title = "Species richness",
                                      facet_experiment = TRUE) {
  .check_flag(facet_experiment, "facet_experiment")
  if (facet_experiment && !is.null(taxonomy_table) && "experiment" %in% names(sequence_table)) {
    if (anyNA(sequence_table$experiment)) stop("experiment must not contain missing values.", call. = FALSE)
    tables <- split(sequence_table, sequence_table$experiment)
    curves <- lapply(names(tables), function(experiment) {
      result <- rarefaction_extrapolation(
        tables[[experiment]], taxonomy_table, n_threads, bin_size, endpoint,
        x_limit, y_limit, show_legend, show_observed, show_extrapolation,
        x_title, y_title, facet_experiment = FALSE
      )
      result$rarefaction_data$experiment <- experiment
      result$rarefaction_data
    })
    curves <- bind_rows(curves)
    colors <- colorRampPalette(c("#3B4992", "#EE0000", "#008B45", "#631879",
      "#008280", "#BB0021", "#5F559B", "#A20056", "#808180", "#1B1919"))(length(unique(curves$tag)))
    names(colors) <- unique(curves$tag)
    fig <- ggplot(curves, aes(depth, richness, colour = tag, group = tag)) +
      geom_line(data = curves[curves$method == "Rarefaction", ], aes(linetype = method), linewidth = 0.7) +
      facet_wrap(~experiment, scales = "free_x") +
      scale_colour_manual(values = colors) +
      scale_linetype_manual(values = c(Rarefaction = "solid", Extrapolation = "dashed")) +
      labs(x = x_title, y = y_title, colour = "Reaction", linetype = NULL) +
      theme_classic(base_size = 10) +
      theme(strip.background = element_blank(),
            legend.position = if (is.null(show_legend)) "none" else show_legend)
    if (show_extrapolation) fig <- fig + geom_line(
      data = curves[curves$method == "Extrapolation", ], aes(linetype = method), linewidth = 0.7)
    if (show_observed) fig <- fig + geom_point(data = curves[curves$method == "Observed", ], size = 1.8)
    if (!is.null(x_limit) || !is.null(y_limit)) fig <- fig + coord_cartesian(xlim = x_limit, ylim = y_limit)
    return(list(rarefaction_data = curves, fig = fig))
  }
  .check_number(n_threads, "n_threads", 1, .Machine$integer.max, TRUE)
  if (!is.null(bin_size)) .check_number(bin_size, "bin_size", 1, Inf, TRUE)
  if (!is.null(endpoint)) .check_number(endpoint, "endpoint", 1, Inf, TRUE)
  .check_flag(show_observed, "show_observed")
  .check_flag(show_extrapolation, "show_extrapolation")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  if (!is.null(taxonomy_table)) {
    .check_columns(sequence_table, "SEQ", "sequence_table")
    .check_columns(taxonomy_table, c("SEQ", "scientific_name"), "taxonomy_table")
    taxonomy_table <- unique(taxonomy_table[c("SEQ", "scientific_name")])
    if (anyDuplicated(taxonomy_table$SEQ)) stop("Each SEQ must have one scientific_name.", call. = FALSE)
    count_columns <- grep("^count_", names(sequence_table), value = TRUE)
    .check_counts(sequence_table[count_columns])
    scientific_name <- taxonomy_table$scientific_name[
      match(sequence_table$SEQ, taxonomy_table$SEQ)
    ]
    keep <- !is.na(scientific_name) & nzchar(scientific_name)
    abundance <- rowsum(
      as.matrix(sequence_table[keep, count_columns, drop = FALSE]),
      scientific_name[keep], reorder = FALSE
    )
  } else {
    abundance <- as.matrix(sequence_table)
  }

  if (!is.numeric(abundance) || any(!is.finite(abundance)) || any(abundance < 0) || any(abundance != floor(abundance))) {
    stop("Rarefaction requires finite non-negative integer read counts.", call. = FALSE)
  }
  abundance <- abundance[, colSums(abundance) > 0, drop = FALSE]
  if (!nrow(abundance) || !ncol(abundance)) stop("No positive read counts remain for rarefaction.", call. = FALSE)
  if (is.null(colnames(abundance))) colnames(abundance) <- paste0("reaction_", seq_len(ncol(abundance)))
  maximum_reads <- max(colSums(abundance))
  if (is.null(endpoint)) endpoint <- ceiling(maximum_reads + maximum_reads * 0.1)

  rarefy_one <- function(counts, tag) {
    counts <- counts[counts > 0]
    total_reads <- sum(counts)
    reaction_bin <- if (is.null(bin_size)) {
      max(1, ceiling(total_reads / 100))
    } else {
      bin_size
    }
    depth <- if (is.null(bin_size)) {
      unique(round(exp(seq(log(1), log(total_reads), length.out = 100))))
    } else {
      unique(c(
        1,
        if (reaction_bin <= total_reads) {
          seq(reaction_bin, total_reads, by = reaction_bin)
        },
        total_reads
      ))
    }
    rarefied_richness <- vapply(depth, function(m) {
      sum(1 - exp(lchoose(total_reads - counts, m) - lchoose(total_reads, m)))
    }, numeric(1))
    curve <- data.frame(
      tag = sub("^count_", "", tag),
      depth = depth,
      richness = rarefied_richness,
      method = "Rarefaction"
    )

    observed <- data.frame(
      tag = sub("^count_", "", tag), depth = total_reads,
      richness = length(counts), method = "Observed"
    )

    if (show_extrapolation) {
      extrapolation_bin <- if (is.null(bin_size)) {
        max(1, ceiling((endpoint - total_reads) / 100))
      } else {
        bin_size
      }
      extrapolation_depth <- unique(c(
        total_reads,
        if (total_reads + extrapolation_bin <= endpoint) {
          seq(total_reads + extrapolation_bin, endpoint, by = extrapolation_bin)
        },
        if (endpoint > total_reads) endpoint
      ))
      f1 <- sum(counts == 1)
      f2 <- sum(counts == 2)
      unseen <- if (f2 > 0) {
        (total_reads - 1) / total_reads * f1^2 / (2 * f2)
      } else {
        (total_reads - 1) / total_reads * f1 * (f1 - 1) / 2
      }
      detection_rate <- if (unseen > 0) f1 / (total_reads * unseen) else 0
      extrapolated_richness <- length(counts) + unseen * (
        1 - (1 - detection_rate)^(extrapolation_depth - total_reads)
      )
      curve <- rbind(
        curve,
        data.frame(
          tag = sub("^count_", "", tag),
          depth = extrapolation_depth,
          richness = extrapolated_richness,
          method = "Extrapolation"
        )
      )
    }
    rbind(curve, observed)
  }

  reaction_index <- seq_len(ncol(abundance))
  calculate_reaction <- function(i) {
    rarefy_one(abundance[, i], colnames(abundance)[i])
  }
  if (n_threads > 1) {
    old_plan <- plan()
    on.exit(plan(old_plan), add = TRUE)
    plan(multisession, workers = n_threads)
    curve_list <- future_lapply(reaction_index, calculate_reaction)
  } else {
    curve_list <- lapply(reaction_index, calculate_reaction)
  }
  rarefaction_data <- bind_rows(curve_list)

  aaas_colors <- c(
    "#3B4992", "#EE0000", "#008B45", "#631879", "#008280",
    "#BB0021", "#5F559B", "#A20056", "#808180", "#1B1919"
  )
  tag_colors <- colorRampPalette(aaas_colors)(ncol(abundance))
  names(tag_colors) <- sub("^count_", "", colnames(abundance))

  plt <- ggplot(mapping = aes(x = depth, y = richness, colour = tag, group = tag) ) +
    geom_line(data = rarefaction_data[rarefaction_data$method == "Rarefaction",,drop = FALSE], aes(linetype = method), alpha = 0.7, linewidth = 0.7) +
    {if(show_extrapolation)
    geom_line(data = rarefaction_data[rarefaction_data$method == "Extrapolation",,drop = FALSE], aes(linetype = method), alpha = 0.7, linewidth = 0.7)
    }+
    {if (show_observed)
      geom_point(
      data = rarefaction_data[rarefaction_data$method == "Observed", , drop = FALSE],
      shape = 19, stroke = 0.5, size = 1.8  )
    } +
    scale_colour_manual(values = tag_colors) +
    scale_linetype_manual(values = c(Rarefaction = "solid", Extrapolation = "dashed")) +
    {if (!is.null(x_limit)) coord_cartesian(xlim = x_limit, ylim = y_limit)} +
    {if (is.null(x_limit) && !is.null(y_limit)) coord_cartesian(ylim = y_limit)} +
    labs(x = x_title, y = y_title, colour = "Reaction", linetype = NULL) +
    theme_classic(base_family = "sans", base_size = 10) +
    theme(
      legend.position = if (is.null(show_legend)) "none" else show_legend,
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.6)
    )

  list(rarefaction_data = rarefaction_data, fig = plt)
}
