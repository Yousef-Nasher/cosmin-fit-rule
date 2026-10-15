#' Re-rate a data set of CFA results under several rules
#'
#' Applies [rate_cfa()] with several rules to each row of a data.frame.
#'
#' @param data A data.frame with one row per CFA result.
#' @param rules Rule names from [cosmin_rules()].
#' @param missing `"reported"` (default) or `"unknown"`; see [rate_cfa()].
#' @param cols Named list giving the column names for `cfi`, `tli`,
#'   `rmsea`, `srmr`, `p` and `df`. Partial lists are filled in with the
#'   defaults. A column that is absent from `data` counts as not reported.
#' @return `data` with one column `rating_<rule>` per rule, plus
#'   `n_indices` (how many of CFI, TLI, RMSEA and SRMR were reported) and
#'   `rmsea_rescue` (TRUE when RMSEA < .06 is the only one of the four
#'   COSMIN clauses that passes among reported indices; NA when RMSEA was
#'   not reported).
#' @examples
#' d <- data.frame(cfi = c(0.97, 0.91, NA), tli = c(NA, 0.89, NA),
#'                 rmsea = c(0.05, 0.05, NA), srmr = c(0.04, 0.09, 0.07),
#'                 p = c(10, 6, 8), df = c(35, 9, 20))
#' rerate(d)
#' rerate(d, rules = c("cosmin", "hu_bentler"), missing = "unknown")
#' @export
rerate <- function(data,
                   rules = c("cosmin", "no_rmsea", "srmr_df", "two_of_three",
                             "all_three"),
                   missing = "reported",
                   cols = list(cfi = "cfi", tli = "tli", rmsea = "rmsea",
                               srmr = "srmr", p = "p", df = "df")) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  missing <- match.arg(missing, c("reported", "unknown"))
  for (r in rules) check_rule(r)
  defaults <- list(cfi = "cfi", tli = "tli", rmsea = "rmsea", srmr = "srmr",
                   p = "p", df = "df")
  cols <- utils::modifyList(defaults, as.list(cols))

  get_col <- function(name) {
    col <- cols[[name]]
    if (!is.null(col) && col %in% names(data)) as_index(data[[col]])
    else rep(NA_real_, nrow(data))
  }
  x <- data.frame(lapply(stats::setNames(names(defaults), names(defaults)),
                         get_col))
  check_ranges(x)
  ind <- pass_indicators(x)

  for (r in rules) {
    data[[paste0("rating_", r)]] <- rate_core(ind, r, missing)
  }
  data$n_indices <- rowSums(!is.na(x[c("cfi", "tli", "rmsea", "srmr")]))
  others <- cbind(ind$cfi_gt95, ind$tli_gt95, ind$srmr_lt08)
  other_pass <- rowSums(!is.na(others) & others) > 0
  data$rmsea_rescue <- ifelse(is.na(ind$rmsea_lt06), NA,
                              ind$rmsea_lt06 & !other_pass)
  data
}
