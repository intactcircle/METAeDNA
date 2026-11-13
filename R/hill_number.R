#' Calculate Hill diversity numbers
#'
#' Calculate Hill numbers for each row or column of an abundance matrix.
#' Common orders are species richness (`q = 0`), the exponential of Shannon
#' entropy (`q = 1`), and inverse Simpson concentration (`q = 2`).
#'
#' @param matrix A numeric abundance matrix or data frame.
#' @param q Non-negative Hill-number order. Common values are 0, 1, and 2.
#' @param margin Use `1` to calculate each row or `2` to calculate each column.
#' @param relative_abundance If `FALSE` (default), each row or column is
#'   converted to relative abundance before calculation. Set to `TRUE` when the
#'   input already contains relative abundances.
#'
#' @return A named numeric vector containing one Hill number per row or column.
#'
#' @examples
#' abundance <- matrix(
#'   c(10, 10, 0, 5, 3, 2),
#'   nrow = 2,
#'   byrow = TRUE,
#'   dimnames = list(c("sample_1", "sample_2"), c("sp_1", "sp_2", "sp_3"))
#' )
#' hill_number(abundance, q = 0)
#' hill_number(abundance, q = 1)
#' hill_number(abundance, q = 2)
#' @export
hill_number <- function(matrix, q = 0, margin = 1,
                        relative_abundance = FALSE) {
  if (!margin %in% c(1, 2)) {
    stop("margin must be 1 for rows or 2 for columns.", call. = FALSE)
  }
  if (length(q) != 1L || !is.numeric(q) || !is.finite(q) || q < 0) {
    stop("q must be one non-negative number.", call. = FALSE)
  }

  matrix <- as.matrix(matrix)
  if (!is.numeric(matrix) || anyNA(matrix) || any(matrix < 0)) {
    stop("matrix must contain non-negative numbers without missing values.",
         call. = FALSE)
  }

  hill_single <- function(abundances) {
    abundances <- abundances[abundances > 0]
    if (!length(abundances)) return(0)
    if (!relative_abundance) abundances <- abundances / sum(abundances)

    if (q == 0) {
      length(abundances)
    } else if (q == 1) {
      exp(-sum(abundances * log(abundances)))
    } else {
      sum(abundances^q)^(1 / (1 - q))
    }
  }

  apply(matrix, margin, hill_single)
}
