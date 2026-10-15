# Generate cross-language test vectors for the Python port of cosminfit.
#
# Usage (from anywhere):  Rscript crosscheck/make_vectors.R
#
# The R package is installed from source (../cosminfit) into a temporary
# library, so the vectors always reflect the current R sources. If that
# install fails, an already installed copy is used instead.
#
# Output (next to this script):
#   test_vectors.csv     rate_cfa() ratings: (a) 81-row pass/fail/NA grid,
#                        (b) 3000 random cases, (c) threshold edges,
#                        (d) string inputs (input cleaning); x 8 rules x 2 modes
#   rerate_vectors.csv   rerate() output (all 8 rules, n_indices,
#                        rmsea_rescue) for sets a-c, both modes
#   mapping_vectors.csv  cfa_df, n_moments, srmr_to_rmsea, rmsea_to_srmr,
#                        rmsea_equivalent_table, wilson_ci, type-7 quantiles,
#                        cluster_boot_prop (deterministic parts)
#
# Numbers are written with sprintf("%.17g") so that every double round-trips
# exactly; NA is written as NA.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
pkg_src <- normalizePath(file.path(here, "..", "cosminfit"), mustWork = TRUE)

lib <- file.path(tempdir(), "cosminfit_lib")
dir.create(lib, showWarnings = FALSE)
status <- system2(file.path(R.home("bin"), "R"),
                  c("CMD", "INSTALL", "--no-test-load", "-l", shQuote(lib),
                    shQuote(pkg_src)), stdout = FALSE, stderr = FALSE)
if (status == 0) {
  library(cosminfit, lib.loc = lib)
  cat("Installed cosminfit", as.character(packageVersion("cosminfit", lib.loc = lib)),
      "from", pkg_src, "into a temporary library\n")
} else {
  message("R CMD INSTALL failed; falling back to an installed cosminfit")
  library(cosminfit)
}

fmt <- function(v) {
  out <- sprintf("%.17g", as.numeric(v))
  out[is.na(v)] <- "NA"
  out
}
rules <- cosmin_rules()$rule
# df formula without cfa_df()'s checks, only to choose identified (p, m)
df_formula <- function(p, m) p * (p + 1) / 2 - (2 * p + m * (m - 1) / 2)
modes <- c("reported", "unknown")

# ---- input sets ---------------------------------------------------------------

# (a) full 3^4 grid: each index passes every variant of its clause (P),
# fails every variant (F), or is not reported (NA). p = 10, df = cfa_df(10).
pf <- list(cfi = c(0.97, 0.90), tli = c(0.97, 0.90),
           rmsea = c(0.04, 0.08), srmr = c(0.03, 0.10))
key <- c("P", "F", NA)
g <- expand.grid(c = key, t = key, r = key, s = key, stringsAsFactors = FALSE)
pick <- function(v, k) ifelse(is.na(k), NA, ifelse(k == "P", v[1], v[2]))
set_a <- data.frame(cfi = pick(pf$cfi, g$c), tli = pick(pf$tli, g$t),
                    rmsea = pick(pf$rmsea, g$r), srmr = pick(pf$srmr, g$s),
                    p = 10, df = cfa_df(10))

# (b) 3000 random cases. Half of the cases are rounded to 3 decimals as in
# published papers (so values hit the cut-offs exactly), half keep full
# double precision. Each index is NA with probability .25, independently.
set.seed(2026)
n <- 3000
p_b <- sample(4:20, n, replace = TRUE)
set_b <- data.frame(cfi = runif(n, 0.85, 1), tli = runif(n, 0.85, 1),
                    rmsea = runif(n, 0, 0.15), srmr = runif(n, 0, 0.15),
                    p = p_b, df = cfa_df(p_b))
rounded <- runif(n) < 0.5
for (v in c("cfi", "tli", "rmsea", "srmr")) {
  set_b[[v]][rounded] <- round(set_b[[v]][rounded], 3)
}
for (v in c("cfi", "tli", "rmsea", "srmr")) {
  set_b[[v]][runif(n) < 0.25] <- NA
}

# (c) values exactly at and next to the thresholds .95, .96, .06, .08, .09.
cut10 <- 0.06 * sqrt(35 / 55)  # srmr_df cut-off for p = 10, df = 35
c1 <- expand.grid(cfi = c(0.95, 0.96, 0.9499, 0.9501, NA),
                  tli = c(0.95, 0.96, NA),
                  rmsea = c(0.06, 0.0599, 0.0601, 0.09, NA),
                  srmr = c(0.08, 0.09, 0.0799, 0.0801, 0.0901, 0.06, cut10, NA))
c1$p <- 10
c1$df <- cfa_df(10)
# p / df problems for the srmr_df rule
pdf <- data.frame(p = c(10, NA, 10, 3, 10, 0, 7.5, NA, 1, 10),
                  df = c(35, 35, NA, 0, -1, 5, 14, NA, 0.5, Inf))
c2 <- merge(expand.grid(cfi = c(0.90, 0.96, NA), tli = NA, rmsea = NA,
                        srmr = c(0.03, 0.05, 0.09, NA)), pdf)
# SRMR exactly at the srmr_df cut-off, for many p and m
pm <- expand.grid(p = 4:30, m = 1:3)
pm <- pm[df_formula(pm$p, pm$m) > 0, ]
pm$df <- cfa_df(pm$p, pm$m)
cutv <- 0.06 * sqrt(pm$df / (pm$p * (pm$p + 1) / 2))
c3 <- data.frame(cfi = NA, tli = NA, rmsea = NA,
                 srmr = c(cutv, cutv * (1 - 1e-15), cutv * (1 + 1e-15)),
                 p = pm$p, df = pm$df)
set_c <- rbind(c1, c2[names(c1)], c3)

# (d) string inputs: input cleaning (NA / non-finite / non-numeric strings
# are "not reported", numeric strings are parsed like as.numeric()).
strs <- c("0.97", " 0.97", "0.97 ", "\t0.97", "n.r.", "", " ", "NA", "NaN",
          "Inf", "-Inf", "inf", "infinity", "1e-1", "9.7E-1", ".97", "+0.97",
          "0x1", "0x1A", "0X1a", "0x1.8p1", "0x", "1_0", "0.97a", "1,0",
          "TRUE", "1e", "1e+", "1.", "1d5", "- 1", "1e5000", "1e-400",
          "0.05", "5e-2", "0.0799", "0.95", "0.06", "0.08", "0.09", "-0.01",
          "1.3", NA)
na_chr <- rep(NA_character_, length(strs))
set_d <- rbind(
  data.frame(cfi = strs, tli = na_chr, rmsea = na_chr, srmr = na_chr),
  data.frame(cfi = na_chr, tli = na_chr, rmsea = strs, srmr = na_chr),
  data.frame(cfi = na_chr, tli = na_chr, rmsea = na_chr, srmr = strs),
  data.frame(cfi = "0.90", tli = strs, rmsea = "0.07", srmr = strs),
  stringsAsFactors = FALSE)
set_d$p <- 10
set_d$df <- 35

sets <- list(a = set_a, b = set_b, c = set_c, d = set_d)

# ---- rate_cfa vectors -----------------------------------------------------------

rows <- list()
for (s in names(sets)) {
  X <- sets[[s]]
  inputs <- data.frame(
    set = s, case = seq_len(nrow(X)),
    cfi = if (s == "d") X$cfi else fmt(X$cfi),
    tli = if (s == "d") X$tli else fmt(X$tli),
    rmsea = if (s == "d") X$rmsea else fmt(X$rmsea),
    srmr = if (s == "d") X$srmr else fmt(X$srmr),
    p = fmt(X$p), df = fmt(X$df), stringsAsFactors = FALSE)
  for (rule in rules) for (mode in modes) {
    rating <- suppressWarnings(
      rate_cfa(X$cfi, X$tli, X$rmsea, X$srmr, rule = rule, p = X$p,
               df = X$df, missing = mode))
    rows[[length(rows) + 1]] <- cbind(inputs, rule = rule, missing = mode,
                                      rating = rating, stringsAsFactors = FALSE)
  }
}
tv <- do.call(rbind, rows)
write.csv(tv, file.path(here, "test_vectors.csv"), row.names = FALSE)
cat("test_vectors.csv:", nrow(tv), "rows;",
    paste(sprintf("%s=%d", names(sets), sapply(sets, nrow)), collapse = " "),
    "cases\n")
print(table(tv$rating))

# ---- rerate vectors -----------------------------------------------------------------

rr <- list()
for (s in c("a", "b", "c")) {
  X <- sets[[s]]
  for (mode in modes) {
    out <- suppressWarnings(rerate(X, rules = rules, missing = mode))
    rr[[length(rr) + 1]] <- data.frame(
      set = s, case = seq_len(nrow(X)), missing = mode,
      cfi = fmt(X$cfi), tli = fmt(X$tli), rmsea = fmt(X$rmsea),
      srmr = fmt(X$srmr), p = fmt(X$p), df = fmt(X$df),
      out[paste0("rating_", rules)], n_indices = out$n_indices,
      rmsea_rescue = out$rmsea_rescue, stringsAsFactors = FALSE)
  }
}
rv <- do.call(rbind, rr)
write.csv(rv, file.path(here, "rerate_vectors.csv"), row.names = FALSE)
cat("rerate_vectors.csv:", nrow(rv), "rows\n")

# ---- mapping / CI vectors --------------------------------------------------------

mv <- list()
add <- function(func, output, value, convention = NA, p = NA, m = NA, df = NA,
                x = NA, n = NA, conf = NA) {
  mv[[length(mv) + 1]] <<- data.frame(
    func = func, convention = convention, p = fmt(p), m = fmt(m),
    df = fmt(df), x = if (is.character(x)) x else fmt(x),
    n = if (is.character(n)) n else fmt(n), conf = fmt(conf),
    output = output, value = fmt(value), stringsAsFactors = FALSE)
}

grid <- expand.grid(p = 4:30, m = 1:3)
grid <- grid[df_formula(grid$p, grid$m) > 0, ]
grid$df <- cfa_df(grid$p, grid$m)
add("cfa_df", "df", grid$df, p = grid$p, m = grid$m)
for (conv in c("lavaan", "offdiag")) {
  add("n_moments", "D", n_moments(1:30, conv), convention = conv, p = 1:30)
  for (s in c(0.08, 0.05, 0.10, 0.0123456789)) {
    add("srmr_to_rmsea", "rmsea", srmr_to_rmsea(s, grid$p, grid$df, conv),
        convention = conv, p = grid$p, m = grid$m, df = grid$df, x = s)
  }
  for (r in c(0.06, 0.05, 0.08)) {
    add("rmsea_to_srmr", "srmr", rmsea_to_srmr(r, grid$p, grid$df, conv),
        convention = conv, p = grid$p, m = grid$m, df = grid$df, x = r)
  }
  for (mm in 1:3) {
    pp <- c(4, 5, 6, 8, 10, 15, 20, 30)
    pp <- pp[df_formula(pp, mm) > 0]
    tab <- rmsea_equivalent_table(p = pp, m = mm, convention = conv)
    add("rmsea_equivalent_table", "rmsea_equivalent", tab$rmsea_equivalent,
        convention = conv, p = tab$p, m = tab$m, df = tab$df, x = 0.08)
    add("rmsea_equivalent_table", "D", tab$D,
        convention = conv, p = tab$p, m = tab$m, df = tab$df, x = 0.08)
  }
}

wx <- expand.grid(xn = c("0/1", "1/1", "0/10", "5/100", "50/100", "3/7",
                         "0/50", "10/50", "50/50", "1/1000", "999/1000",
                         "17/23", "123/456", "2.5/10"),
                  conf = c(0.95, 0.90, 0.99, 0.80, 0.50, 0.999),
                  stringsAsFactors = FALSE)
wx$x <- as.numeric(sub("/.*", "", wx$xn))
wx$n <- as.numeric(sub(".*/", "", wx$xn))
for (i in seq_len(nrow(wx))) {
  w <- wilson_ci(wx$x[i], wx$n[i], conf = wx$conf[i])
  for (o in c("est", "lower", "upper")) {
    add("wilson_ci", o, w[[o]], x = wx$x[i], n = wx$n[i], conf = wx$conf[i])
  }
  add("qnorm", "z", qnorm(1 - (1 - wx$conf[i]) / 2), conf = wx$conf[i])
}

# R's default (type 7) quantile, as used by cluster_boot_prop()
set.seed(11)
for (len in c(1, 2, 3, 7, 50, 501)) {
  v <- round(runif(len), sample(c(2, 17), 1))
  for (pr in c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1)) {
    add("quantile_type7", "q", quantile(v, pr, names = FALSE),
        x = paste(fmt(v), collapse = ";"), conf = pr)
  }
}

# cluster_boot_prop: deterministic parts (est, n, n_clusters) and an
# R interval with large B for a Monte Carlo comparison.
set.seed(42)
bx <- rbinom(60, 1, 0.4) == 1
bx[c(3, 17)] <- NA
bcl <- rep(1:20, each = 3)
bres <- cluster_boot_prop(bx, bcl, B = 20000, seed = 7)
xs <- paste(ifelse(is.na(bx), "NA", as.integer(bx)), collapse = ";")
cs <- paste(bcl, collapse = ";")
for (o in c("est", "lower", "upper", "B", "n", "n_clusters")) {
  add("cluster_boot_prop", o, bres[[o]], x = xs, n = cs, conf = 0.95)
}

mvd <- do.call(rbind, mv)
write.csv(mvd, file.path(here, "mapping_vectors.csv"), row.names = FALSE)
cat("mapping_vectors.csv:", nrow(mvd), "rows\n")
print(table(mvd$func))
