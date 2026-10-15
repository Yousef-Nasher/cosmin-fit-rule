# Each rule is a small logic tree over "pass indicators" (leaves).
# Nodes: or(...), and(...), at_least(k, ...). Leaves are indicator names.
node_or <- function(...) list(op = "or", args = list(...))
node_and <- function(...) list(op = "and", args = list(...))
node_at_least <- function(k, ...) list(op = "at_least", k = k, args = list(...))

rule_trees <- function() {
  list(
    cosmin = node_or("cfi_gt95", "tli_gt95", "rmsea_lt06", "srmr_lt08"),
    no_rmsea = node_or("cfi_gt95", "tli_gt95", "srmr_lt08"),
    srmr_df = node_or("cfi_gt95", "tli_gt95", "rmsea_lt06", "srmr_lt_df"),
    two_of_three = node_at_least(2, "cfi_gt95", "rmsea_lt06", "srmr_lt08"),
    all_three = node_and("cfi_gt95", "rmsea_lt06", "srmr_lt08"),
    hu_bentler = node_or(node_and("cfi_ge95", "srmr_le09"),
                         node_and("rmsea_le06", "srmr_le09")),
    prinsen2016_a = node_and(node_or("cfi_gt95", "tli_gt95"),
                             node_or("rmsea_lt06", "srmr_lt08")),
    prinsen2016_b = node_or(node_and(node_or("cfi_gt95", "tli_gt95"),
                                     "rmsea_lt06"),
                            "srmr_lt08")
  )
}

#' Available rating rules
#'
#' Lists the rules that [rate_cfa()] and [rerate()] can apply. All
#' inequalities are as written (strict unless shown otherwise).
#'
#' * `cosmin`: CFI > .95 OR TLI > .95 OR RMSEA < .06 OR SRMR < .08
#'   (Prinsen et al. 2018; Mokkink et al. 2024).
#' * `no_rmsea`: CFI > .95 OR TLI > .95 OR SRMR < .08.
#' * `srmr_df`: CFI > .95 OR TLI > .95 OR RMSEA < .06 OR
#'   SRMR < .06 * sqrt(df / D), with D = p(p+1)/2. Needs `p` and `df`.
#' * `two_of_three`: at least two of CFI > .95, RMSEA < .06, SRMR < .08
#'   (TLI not used).
#' * `all_three`: CFI > .95 AND RMSEA < .06 AND SRMR < .08.
#' * `hu_bentler`: (CFI >= .95 AND SRMR <= .09) OR
#'   (RMSEA <= .06 AND SRMR <= .09).
#' * `prinsen2016_a`, `prinsen2016_b`: two readings of the 2016
#'   COSMIN-COMET guideline (Prinsen et al. 2016, Trials 17:449, Table 4),
#'   which printed "CFI or TLI or comparable measure > 0.95 AND RMSEA < 0.06
#'   OR SRMR < 0.08". Reading a is (CFI > .95 OR TLI > .95) AND
#'   (RMSEA < .06 OR SRMR < .08); reading b is
#'   ((CFI > .95 OR TLI > .95) AND RMSEA < .06) OR SRMR < .08.
#'
#' @return A data.frame with columns `rule`, `definition`, `indices` (the
#'   indices the rule uses) and `needs_p_df` (whether `p` and `df` are
#'   needed).
#' @examples
#' cosmin_rules()
#' @export
cosmin_rules <- function() {
  data.frame(
    rule = names(rule_trees()),
    definition = c(
      "CFI > .95 OR TLI > .95 OR RMSEA < .06 OR SRMR < .08",
      "CFI > .95 OR TLI > .95 OR SRMR < .08",
      "CFI > .95 OR TLI > .95 OR RMSEA < .06 OR SRMR < .06*sqrt(df/D)",
      "at least two of: CFI > .95, RMSEA < .06, SRMR < .08",
      "CFI > .95 AND RMSEA < .06 AND SRMR < .08",
      "(CFI >= .95 AND SRMR <= .09) OR (RMSEA <= .06 AND SRMR <= .09)",
      "(CFI > .95 OR TLI > .95) AND (RMSEA < .06 OR SRMR < .08)",
      "((CFI > .95 OR TLI > .95) AND RMSEA < .06) OR SRMR < .08"
    ),
    indices = c("cfi, tli, rmsea, srmr", "cfi, tli, srmr",
                "cfi, tli, rmsea, srmr", "cfi, rmsea, srmr",
                "cfi, rmsea, srmr", "cfi, rmsea, srmr",
                "cfi, tli, rmsea, srmr", "cfi, tli, rmsea, srmr"),
    needs_p_df = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

#' Rate CFA fit against COSMIN and alternative rules
#'
#' Rates each set of fit indices as `"sufficient"`, `"insufficient"` or
#' `"indeterminate"` under one rule (see [cosmin_rules()]). Inputs are
#' recycled like data.frame columns. Values that are `NA`, not finite or
#' cannot be read as a number count as "not reported".
#'
#' Each rule is a logic formula over pass indicators (for example
#' `CFI > .95`). An indicator is missing when its index is not reported. For
#' `srmr_df`, the SRMR term is missing also when `p` or `df` is missing or
#' `df <= 0`.
#'
#' **`missing = "reported"`** (primary analysis: use only the indices the
#' study reported). The formula is evaluated bottom-up; each part is either
#' evaluable (TRUE/FALSE) or not:
#' * an OR group is evaluable if at least one member is evaluable; it is
#'   TRUE if any evaluable member is TRUE;
#' * an AND group is evaluable only if all members are evaluable;
#' * "at least k of" is evaluable if at least k members are evaluable; it
#'   is TRUE if at least k evaluable members are TRUE.
#'
#' The rating is sufficient (TRUE), insufficient (FALSE) or indeterminate
#' (not evaluable). So for `cosmin`, `no_rmsea` and `srmr_df` any reported
#' passing index gives sufficient; `two_of_three` needs two reported indices;
#' `all_three` needs all three; each `hu_bentler` conjunct needs both its
#' indices. `prinsen2016_a` needs at least one of CFI/TLI and at least one of
#' RMSEA/SRMR. `prinsen2016_b` is evaluable if SRMR is reported, or if RMSEA
#' and at least one of CFI/TLI are reported. For example, CFI = .97 and
#' SRMR = .10 with RMSEA not reported is insufficient: the part
#' "(CFI or TLI) and RMSEA" cannot be evaluated and is ignored, and the SRMR
#' part fails. Under `missing = "unknown"` the same input is indeterminate.
#'
#' **`missing = "unknown"`** (sensitivity analysis). A missing indicator
#' may be TRUE or FALSE. All completions are enumerated (at most 2^4); the
#' rating is sufficient or insufficient only if every completion gives the
#' same result, otherwise indeterminate. This equals three-valued (Kleene)
#' logic.
#'
#' With all indices of a rule reported, both modes give the same rating.
#'
#' @param cfi,tli,rmsea,srmr Fit index values.
#' @param rule One rule name from [cosmin_rules()]. Default `"cosmin"`.
#' @param p Number of items (only used by `srmr_df`).
#' @param df Model degrees of freedom (only used by `srmr_df`).
#' @param missing How to treat indices that were not reported:
#'   `"reported"` (default) or `"unknown"`.
#' @return Character vector of ratings.
#' @examples
#' rate_cfa(cfi = 0.93, rmsea = 0.05)
#' rate_cfa(cfi = 0.93, rmsea = 0.07, srmr = 0.05, rule = "two_of_three")
#' rate_cfa(cfi = 0.97, rule = "all_three")
#' rate_cfa(cfi = 0.93, rule = "cosmin", missing = "unknown")
#'
#' # Several studies at once
#' rate_cfa(cfi = c(0.96, 0.90, NA), srmr = c(NA, 0.09, NA))
#'
#' # From a lavaan fit
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   fit <- lavaan::cfa("visual =~ x1 + x2 + x3 + x4 + x5 + x6",
#'                      data = lavaan::HolzingerSwineford1939)
#'   fm <- lavaan::fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr", "df"))
#'   rate_cfa(cfi = fm["cfi"], tli = fm["tli"], rmsea = fm["rmsea"],
#'            srmr = fm["srmr"], rule = "srmr_df", p = 6, df = fm["df"])
#' }
#' @export
rate_cfa <- function(cfi = NA, tli = NA, rmsea = NA, srmr = NA,
                     rule = "cosmin", p = NA, df = NA,
                     missing = c("reported", "unknown")) {
  missing <- match.arg(missing)
  rule <- check_rule(rule)
  args <- list(cfi = cfi, tli = tli, rmsea = rmsea, srmr = srmr, p = p, df = df)
  if (any(lengths(args) == 0)) return(character(0))
  x <- as.data.frame(lapply(args, as_index))
  check_ranges(x)
  rate_core(pass_indicators(x), rule, missing)
}

# ---- internal engine -------------------------------------------------------

check_rule <- function(rule) {
  ok <- is.character(rule) && length(rule) == 1 && rule %in% names(rule_trees())
  if (!ok) {
    stop("`rule` must be one of: ",
         paste(names(rule_trees()), collapse = ", "), call. = FALSE)
  }
  rule
}

# Coerce to numeric; anything not a finite number becomes NA.
as_index <- function(x) {
  x <- unname(x)
  if (is.factor(x)) x <- as.character(x)
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA
  x
}

check_ranges <- function(x) {
  if (any(x$cfi > 1.2, x$tli > 1.2, na.rm = TRUE)) {
    warning("Some CFI or TLI values are above 1.2; check the input.",
            call. = FALSE)
  }
  if (any(x$rmsea < 0, x$srmr < 0, na.rm = TRUE)) {
    warning("Some RMSEA or SRMR values are negative; check the input.",
            call. = FALSE)
  }
  invisible(TRUE)
}

# Logical pass indicators; NA means "not reported".
pass_indicators <- function(x) {
  ok_pdf <- !is.na(x$p) & !is.na(x$df) & x$p >= 1 & x$df > 0
  srmr_cut_df <- rep(NA_real_, nrow(x))
  srmr_cut_df[ok_pdf] <- 0.06 * sqrt(x$df[ok_pdf] / (x$p[ok_pdf] * (x$p[ok_pdf] + 1) / 2))
  list(
    cfi_gt95 = x$cfi > 0.95,
    tli_gt95 = x$tli > 0.95,
    rmsea_lt06 = x$rmsea < 0.06,
    srmr_lt08 = x$srmr < 0.08,
    srmr_lt_df = x$srmr < srmr_cut_df,
    cfi_ge95 = x$cfi >= 0.95,
    rmsea_le06 = x$rmsea <= 0.06,
    srmr_le09 = x$srmr <= 0.09
  )
}

rate_core <- function(ind, rule, missing) {
  tree <- rule_trees()[[rule]]
  res <- if (missing == "reported") eval_reported(tree, ind) else eval_unknown(tree, ind)
  out <- rep("indeterminate", length(res))
  out[res %in% TRUE] <- "sufficient"
  out[res %in% FALSE] <- "insufficient"
  out
}

# Evaluate a tree with "reported" semantics. NA = not evaluable.
# With no NA in the leaves this is ordinary Boolean logic.
eval_reported <- function(node, ind) {
  if (is.character(node)) return(ind[[node]])
  vals <- do.call(cbind, lapply(node$args, eval_reported, ind = ind))
  n_eval <- rowSums(!is.na(vals))
  n_true <- rowSums(!is.na(vals) & vals)
  res <- switch(node$op,
    or = n_true > 0,
    and = n_true == ncol(vals),
    at_least = n_true >= node$k
  )
  evaluable <- switch(node$op,
    or = n_eval > 0,
    and = n_eval == ncol(vals),
    at_least = n_eval >= node$k
  )
  res[!evaluable] <- NA
  res
}

# Evaluate a tree with "unknown" semantics: try every TRUE/FALSE value of the
# missing leaves; keep the result only if it is the same for all of them.
eval_unknown <- function(node, ind) {
  leaves <- unique(tree_leaves(node))
  grid <- expand.grid(rep(list(c(FALSE, TRUE)), length(leaves)))
  results <- sapply(seq_len(nrow(grid)), function(i) {
    filled <- ind
    for (j in seq_along(leaves)) {
      v <- filled[[leaves[j]]]
      v[is.na(v)] <- grid[i, j]
      filled[[leaves[j]]] <- v
    }
    eval_reported(node, filled)
  })
  results <- matrix(results, ncol = nrow(grid))
  res <- rep(NA, nrow(results))
  res[rowSums(results) == ncol(results)] <- TRUE
  res[rowSums(results) == 0] <- FALSE
  res
}

tree_leaves <- function(node) {
  if (is.character(node)) return(node)
  unlist(lapply(node$args, tree_leaves))
}
