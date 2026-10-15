#' Degrees of freedom of a simple-structure CFA
#'
#' Degrees of freedom of a CFA model with `p` items loading on `m`
#' correlated factors, each item on exactly one factor:
#' `p * (p + 1) / 2 - (2 * p + m * (m - 1) / 2)`.
#'
#' @param p Number of items (vector allowed).
#' @param m Number of factors (vector allowed, recycled with `p`).
#' @return Integer-valued numeric vector of degrees of freedom. An error is
#'   raised if any df is negative. A df of 0 gives a warning, because the
#'   SRMR-RMSEA mapping is then undefined.
#' @examples
#' cfa_df(4)
#' cfa_df(10, m = 2)
#' cfa_df(c(4, 6, 8, 10))
#' @export
cfa_df <- function(p, m = 1) {
  check_whole(p, "p")
  check_whole(m, "m")
  df <- p * (p + 1) / 2 - (2 * p + m * (m - 1) / 2)
  if (any(df < 0)) {
    stop("The model is not identified (negative df) for some p and m.",
         call. = FALSE)
  }
  if (any(df == 0)) {
    warning("df = 0 for some p and m; the SRMR-RMSEA mapping is undefined.",
            call. = FALSE)
  }
  df
}

#' Number of sample moments used by SRMR
#'
#' The number of elements `D` averaged by SRMR for `p` items.
#'
#' * `"lavaan"`: `D = p (p + 1) / 2`, including the diagonal. This is
#'   the convention of lavaan and of Shi, Maydeu-Olivares and Rosseel (2020),
#'   who divide by \eqn{p(p+1)/2} even for correlation structures.
#' * `"offdiag"`: `D = p (p - 1) / 2`, correlations below the diagonal only
#'   (an alternative, for sensitivity checks).
#'
#' @param p Number of items (vector allowed).
#' @param convention `"lavaan"` (default) or `"offdiag"`.
#' @return Numeric vector `D`.
#' @examples
#' n_moments(10)
#' n_moments(10, convention = "offdiag")
#' @export
n_moments <- function(p, convention = c("lavaan", "offdiag")) {
  convention <- match.arg(convention)
  check_whole(p, "p")
  if (convention == "lavaan") p * (p + 1) / 2 else p * (p - 1) / 2
}

# Stop unless x is a non-empty vector of positive whole numbers.
check_whole <- function(x, name) {
  ok <- is.numeric(x) && length(x) > 0 && all(is.finite(x)) &&
    all(x >= 1) && all(x == round(x))
  if (!ok) stop("`", name, "` must contain positive whole numbers.",
                call. = FALSE)
  invisible(TRUE)
}
