#' Population relation between SRMR and RMSEA
#'
#' Let `G` be the sum of squared residual correlations. Then
#' `SRMR = sqrt(G / D)` and, at population level under unweighted least
#' squares (ULS), `RMSEA = sqrt(G / df)`. Hence
#' `RMSEA = SRMR * sqrt(D / df)` and `SRMR = RMSEA * sqrt(df / D)`, where
#' `D` is given by [n_moments()].
#'
#' This holds for population values under ULS. It does not hold exactly for
#' sample values or for other estimators (such as ML), so do not use it to
#' convert one reported index into another.
#'
#' @param srmr,rmsea Index values (vectors allowed).
#' @param p Number of items.
#' @param df Model degrees of freedom (must be positive), e.g. from
#'   [cfa_df()].
#' @param convention Passed to [n_moments()]: `"lavaan"` (default) or
#'   `"offdiag"`.
#' @return Numeric vector. `NA` (with a warning) where `df <= 0`.
#' @examples
#' srmr_to_rmsea(0.08, p = 4, df = 2)
#' rmsea_to_srmr(0.06, p = 10, df = cfa_df(10))
#' @export
srmr_to_rmsea <- function(srmr, p, df, convention = c("lavaan", "offdiag")) {
  convention <- match.arg(convention)
  srmr * sqrt(n_moments(p, convention) / positive_df(df))
}

#' @rdname srmr_to_rmsea
#' @export
rmsea_to_srmr <- function(rmsea, p, df, convention = c("lavaan", "offdiag")) {
  convention <- match.arg(convention)
  rmsea * sqrt(positive_df(df) / n_moments(p, convention))
}

#' RMSEA equivalent of an SRMR cut-off
#'
#' For each number of items, the RMSEA value that corresponds, at population
#' level under ULS, to the SRMR cut-off in a one- or multi-factor
#' simple-structure CFA. See [srmr_to_rmsea()].
#'
#' @param p Numbers of items.
#' @param m Number of factors (recycled with `p`).
#' @param srmr_cut SRMR cut-off. Default 0.08 (COSMIN).
#' @param convention Passed to [n_moments()].
#' @return A data.frame with columns `p`, `m`, `df`, `D` and
#'   `rmsea_equivalent`.
#' @examples
#' rmsea_equivalent_table()
#' rmsea_equivalent_table(p = c(6, 9, 12), m = 3)
#' @export
rmsea_equivalent_table <- function(p = c(4, 5, 6, 8, 10, 15, 20, 30), m = 1,
                                   srmr_cut = 0.08,
                                   convention = c("lavaan", "offdiag")) {
  convention <- match.arg(convention)
  out <- data.frame(p = p, m = m)
  out$df <- cfa_df(out$p, out$m)
  out$D <- n_moments(out$p, convention)
  out$rmsea_equivalent <- srmr_to_rmsea(srmr_cut, out$p, out$df, convention)
  out
}

# Replace non-positive df by NA, with a warning.
positive_df <- function(df) {
  bad <- !is.na(df) & df <= 0
  if (any(bad)) {
    warning("df must be positive; returning NA where df <= 0.", call. = FALSE)
    df[bad] <- NA
  }
  df
}
