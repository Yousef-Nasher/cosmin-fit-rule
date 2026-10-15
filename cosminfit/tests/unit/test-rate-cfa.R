# Fixed values: P = passes every variant of its clause, F = fails every variant.
# p = 10, df = 35 gives a scaled SRMR cut-off of .06 * sqrt(35/55) = .0479.
vals <- list(cfi = c(P = 0.97, F = 0.90), tli = c(P = 0.97, F = 0.90),
             rmsea = c(P = 0.04, F = 0.08), srmr = c(P = 0.03, F = 0.10))
rules <- cosmin_rules()$rule
check(identical(rules, c("cosmin", "no_rmsea", "srmr_df", "two_of_three",
                         "all_three", "hu_bentler", "prinsen2016_a",
                         "prinsen2016_b")), "rule names")

# Independent reference, written directly from the specification.
# c, t, r, s are TRUE (reported, pass), FALSE (reported, fail) or NA.
any_rep <- function(...) { v <- c(...); if (all(is.na(v))) NA else any(v, na.rm = TRUE) }
ref_reported <- list(
  cosmin = function(c, t, r, s) any_rep(c, t, r, s),
  no_rmsea = function(c, t, r, s) any_rep(c, t, s),
  srmr_df = function(c, t, r, s) any_rep(c, t, r, s),
  two_of_three = function(c, t, r, s) {
    v <- c(c, r, s); if (sum(!is.na(v)) < 2) NA else sum(v, na.rm = TRUE) >= 2
  },
  all_three = function(c, t, r, s) if (anyNA(c(c, r, s))) NA else c && r && s,
  hu_bentler = function(c, t, r, s) {
    a <- if (anyNA(c(c, s))) NA else c && s
    b <- if (anyNA(c(r, s))) NA else r && s
    any_rep(a, b)
  },
  prinsen2016_a = function(c, t, r, s) {
    g1 <- any_rep(c, t); g2 <- any_rep(r, s)
    if (anyNA(c(g1, g2))) NA else g1 && g2
  },
  prinsen2016_b = function(c, t, r, s) {
    g1 <- any_rep(c, t)
    a <- if (anyNA(c(g1, r))) NA else g1 && r
    any_rep(a, s)
  }
)
ref_unknown <- function(f, c, t, r, s) {
  x <- c(c, t, r, s)
  miss <- which(is.na(x))
  out <- c()
  for (i in seq_len(2^length(miss))) {
    y <- x
    y[miss] <- as.logical(bitwAnd(i - 1, 2^(seq_along(miss) - 1)))
    out <- c(out, f(y[1], y[2], y[3], y[4]))
  }
  if (all(out)) TRUE else if (!any(out)) FALSE else NA
}
lab <- function(v) ifelse(is.na(v), "indeterminate", ifelse(v, "sufficient", "insufficient"))

# Full 3^4 grid: every index pass, fail or not reported.
g <- expand.grid(c = c(TRUE, FALSE, NA), t = c(TRUE, FALSE, NA),
                 r = c(TRUE, FALSE, NA), s = c(TRUE, FALSE, NA))
pick <- function(v, key) ifelse(is.na(key), NA, ifelse(key, v["P"], v["F"]))
num <- data.frame(cfi = pick(vals$cfi, g$c), tli = pick(vals$tli, g$t),
                  rmsea = pick(vals$rmsea, g$r), srmr = pick(vals$srmr, g$s))
for (rule in rules) {
  for (mode in c("reported", "unknown")) {
    got <- rate_cfa(num$cfi, num$tli, num$rmsea, num$srmr, rule = rule,
                    p = 10, df = 35, missing = mode)
    want <- vapply(seq_len(nrow(g)), function(i) {
      f <- ref_reported[[rule]]
      if (mode == "reported") f(g$c[i], g$t[i], g$r[i], g$s[i])
      else ref_unknown(f, g$c[i], g$t[i], g$r[i], g$s[i])
    }, logical(1))
    check(identical(got, lab(want)), paste("81-row truth table:", rule, mode))
  }
}

# All indices missing -> indeterminate for every rule and mode.
for (rule in rules) for (mode in c("reported", "unknown")) {
  check(rate_cfa(rule = rule, missing = mode) == "indeterminate",
        paste("all NA:", rule, mode))
}

# Single index reported.
check(rate_cfa(cfi = 0.97) == "sufficient", "cosmin, CFI pass alone")
check(rate_cfa(cfi = 0.90) == "insufficient", "cosmin, CFI fail alone (reported)")
check(rate_cfa(cfi = 0.90, missing = "unknown") == "indeterminate", "cosmin, CFI fail alone (unknown)")
check(rate_cfa(rmsea = 0.05) == "sufficient", "cosmin, RMSEA pass alone")
check(rate_cfa(rmsea = 0.05, rule = "no_rmsea") == "indeterminate", "no_rmsea ignores RMSEA")
check(rate_cfa(cfi = 0.97, rule = "two_of_three") == "indeterminate", "two_of_three, one index")
check(rate_cfa(tli = 0.97, srmr = 0.03, rule = "two_of_three") == "indeterminate", "two_of_three ignores TLI")
check(rate_cfa(cfi = 0.97, rule = "all_three") == "indeterminate", "all_three, one index (reported)")
check(rate_cfa(cfi = 0.90, rule = "all_three", missing = "unknown") == "insufficient", "all_three, one fail (unknown)")
check(rate_cfa(cfi = 0.97, rmsea = 0.04, rule = "two_of_three", missing = "unknown") == "sufficient", "two_of_three, two pass (unknown)")

# Equality edges.
check(rate_cfa(cfi = 0.95) == "insufficient", "CFI exactly .95 is not > .95")
check(rate_cfa(tli = 0.95) == "insufficient", "TLI exactly .95 is not > .95")
check(rate_cfa(rmsea = 0.06) == "insufficient", "RMSEA exactly .06 fails")
check(rate_cfa(srmr = 0.08) == "insufficient", "SRMR exactly .08 fails for cosmin")
check(rate_cfa(srmr = 0.0799) == "sufficient", "SRMR .0799 passes")
check(rate_cfa(cfi = 0.95, srmr = 0.09, rule = "hu_bentler") == "sufficient", "hu_bentler CFI .95 and SRMR .09")
check(rate_cfa(rmsea = 0.06, srmr = 0.09, rule = "hu_bentler") == "sufficient", "hu_bentler RMSEA .06 and SRMR .09")
check(rate_cfa(cfi = 0.97, srmr = 0.0901, rule = "hu_bentler") == "insufficient", "hu_bentler SRMR above .09")
check(rate_cfa(cfi = 0.97, rule = "hu_bentler") == "indeterminate", "hu_bentler conjunct needs both")
check(rate_cfa(cfi = 0.90, rmsea = 0.04, srmr = 0.09, rule = "hu_bentler") == "sufficient", "hu_bentler second conjunct")

# srmr_df needs p and df for the SRMR term.
check(rate_cfa(srmr = 0.03, rule = "srmr_df") == "indeterminate", "srmr_df without p/df")
check(rate_cfa(srmr = 0.03, rule = "srmr_df", p = 10) == "indeterminate", "srmr_df without df")
check(rate_cfa(srmr = 0.03, rule = "srmr_df", p = 10, df = 35) == "sufficient", "srmr_df pass")
check(rate_cfa(srmr = 0.06, rule = "srmr_df", p = 10, df = 35) == "insufficient", "srmr_df fail (.06 > .0479)")
check(rate_cfa(srmr = 0.06) == "sufficient", "same SRMR passes cosmin")
check(rate_cfa(srmr = 0.03, rule = "srmr_df", p = 3, df = 0) == "indeterminate", "srmr_df with df = 0")
check(rate_cfa(cfi = 0.97, rule = "srmr_df") == "sufficient", "srmr_df CFI pass without p/df")

# Prinsen 2016 readings.
check(rate_cfa(cfi = 0.97, srmr = 0.03, rule = "prinsen2016_a") == "sufficient", "p16a pass")
check(rate_cfa(cfi = 0.97, rule = "prinsen2016_a") == "indeterminate", "p16a needs RMSEA or SRMR")
check(rate_cfa(cfi = 0.90, srmr = 0.03, rule = "prinsen2016_a") == "insufficient", "p16a CFI fails")
check(rate_cfa(srmr = 0.03, rule = "prinsen2016_b") == "sufficient", "p16b SRMR alone passes")
check(rate_cfa(cfi = 0.97, srmr = 0.10, rule = "prinsen2016_b") == "insufficient", "p16b SRMR fail, no RMSEA")
check(rate_cfa(cfi = 0.97, srmr = 0.10, rule = "prinsen2016_b", missing = "unknown") == "indeterminate", "p16b SRMR fail, RMSEA unknown")
check(rate_cfa(cfi = 0.97, rule = "prinsen2016_b") == "indeterminate", "p16b CFI alone")
check(rate_cfa(cfi = 0.97, rmsea = 0.04, srmr = 0.10, rule = "prinsen2016_b") == "sufficient", "p16b CFI and RMSEA")
check(rate_cfa(cfi = 0.97, rmsea = 0.08, rule = "prinsen2016_b") == "insufficient", "p16b RMSEA fails, no SRMR")
check(rate_cfa(cfi = 0.97, rmsea = 0.08, rule = "prinsen2016_b", missing = "unknown") == "indeterminate", "p16b unknown SRMR")

# Input handling.
check(identical(rate_cfa(cfi = c(0.97, 0.90), srmr = 0.10),
                c("sufficient", "insufficient")), "recycling")
check(length(rate_cfa(cfi = numeric(0))) == 0, "zero length input")
check_error(rate_cfa(cfi = c(0.9, 0.9), srmr = c(0.1, 0.1, 0.1)), "incompatible lengths")
check(rate_cfa(cfi = "n.r.", srmr = 0.05) == "sufficient", "non-numeric = not reported")
check(rate_cfa(cfi = "n.r.") == "indeterminate", "non-numeric alone")
check(rate_cfa(cfi = "0.97") == "sufficient", "numeric string is read")
check(rate_cfa(cfi = Inf) == "indeterminate", "Inf = not reported")
check(rate_cfa(cfi = c(a = 0.97)) == "sufficient", "named input")
check_error(rate_cfa(cfi = 0.97, rule = "nope"), "unknown rule")
check_error(rate_cfa(cfi = 0.97, missing = "nope"), "unknown missing mode")
check_warning(rate_cfa(cfi = 1.3), "CFI above 1.2 warns")
check_warning(rate_cfa(tli = 1.25), "TLI above 1.2 warns")
check_warning(rate_cfa(srmr = -0.01), "negative SRMR warns")
check_warning(rate_cfa(rmsea = -0.01), "negative RMSEA warns")
check_no_warning(rate_cfa(cfi = 1.1, tli = 1.15, rmsea = 0, srmr = 0), "valid edges do not warn")
