#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter mutate arrange group_by summarise distinct rename case_when join_by left_join full_join right_join inner_join anti_join
#' @importFrom tidyr pivot_longer pivot_wider separate unite drop_na
#' @importFrom readr write_csv read_csv
#' @import stringr
#' @import jsonlite
#' @import curl
#' @useDynLib METAeDNA, .registration = TRUE
#' @import Rcpp
NULL
