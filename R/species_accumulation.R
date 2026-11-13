#' Calculate species accumulation across sampling sites
#'
#' Calculate expected richness for increasing numbers of randomly selected
#' sites without replacement. Abundance values are converted to presence or
#' absence. Each matrix row must represent one sampling site, not a PCR replicate.
#'
#' @param community_matrix Numeric site-by-taxon matrix. Counts and relative
#'   abundances are accepted; values greater than zero indicate detection.
#' @param experiment Optional experiment labels in the same order as matrix
#'   rows. Curves are calculated independently and displayed in separate panels.
#'   The default `NULL` calculates one curve from all sites.
#' @param taxa_in_rows If `TRUE`, taxa are rows and sites are columns.
#' @return A list with `accumulation_data` and the ggplot object `fig`.
#'   The table contains `experiment`, `sites` (number sampled), `richness`
#'   (expected number of detected taxa), and `additional_richness` (expected
#'   increase from adding one site within the observed sampling range).
#' @details
#' For N sites and a taxon detected at f sites, its probability of detection
#' in m randomly selected sites is `1 - choose(N - f, m) / choose(N, m)`.
#' Expected richness is the sum of these probabilities across taxa. The curve
#' includes empty sites and ends at the observed number of sites. It does not
#' extrapolate to unsampled sites or estimate undetected species.
#'
#' A flattening curve indicates diminishing gains within the sampled sites;
#' it does not prove that the regional species pool has been fully sampled.
#' Experiment labels define independent groups, not comparable sampling units
#' when different markers or survey methods detect different taxon sets.
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs
#'   theme_classic theme element_blank
#' @examples
#' set.seed(12)
#' community <- matrix(rbinom(20 * 30, 1, 0.15), nrow = 20)
#' groups <- rep(c("experiment_1", "experiment_2"), each = 10)
#' result <- species_accumulation(community, experiment = groups)
#' result$fig
#' tail(result$accumulation_data)
#' @export
species_accumulation <- function(community_matrix, experiment = NULL,
                                 taxa_in_rows = FALSE) {
  .check_flag(taxa_in_rows, "taxa_in_rows")
  abundance <- as.matrix(community_matrix)
  if (taxa_in_rows) abundance <- t(abundance)
  if (!is.numeric(abundance) || any(!is.finite(abundance)) || any(abundance < 0) ||
      !nrow(abundance) || !ncol(abundance)) {
    stop("community_matrix must contain finite non-negative abundances and at least one site and taxon.", call. = FALSE)
  }
  if (is.null(experiment)) experiment <- rep("All sites", nrow(abundance))
  if (length(experiment) != nrow(abundance) || anyNA(experiment) || any(!nzchar(as.character(experiment)))) {
    stop("experiment must supply one non-missing label per site.", call. = FALSE)
  }
  groups <- split(seq_len(nrow(abundance)), as.character(experiment))
  curves <- lapply(names(groups), function(group) {
    rows <- groups[[group]]
    n <- length(rows)
    frequency <- colSums(abundance[rows, , drop = FALSE] > 0)
    frequency <- frequency[frequency > 0]
    richness <- vapply(0:n, function(m) {
      sum(-expm1(lchoose(n - frequency, m) - lchoose(n, m)))
    }, numeric(1))
    data.frame(experiment = group, sites = 0:n, richness = richness,
               additional_richness = c(NA_real_, diff(richness)))
  })
  accumulation_data <- do.call(rbind, curves)
  fig <- ggplot(accumulation_data, aes(sites, richness)) +
    geom_line(colour = "#3B4992", linewidth = 0.7) +
    geom_point(data = accumulation_data[accumulation_data$sites > 0, ], colour = "#3B4992", size = 1.5) +
    facet_wrap(~experiment, scales = "free_x") +
    labs(x = "Number of sites", y = "Expected taxon richness") +
    theme_classic(base_size = 10) + theme(strip.background = element_blank())
  list(accumulation_data = accumulation_data, fig = fig)
}
