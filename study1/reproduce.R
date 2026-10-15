# reproduce.R -- Study 1: re-run the population analysis and compare it with the stored
# results. Run from study1/:
#   Rscript reproduce.R           part A (ULS identities and Table 1; a few seconds)
#   Rscript reproduce.R --dwls    also part B (48 samples of N = 200,000; about 15 minutes)
# Needs lavaan and the local package cosminfit. Values are compared with a numerical
# tolerance, because lavaan's optimiser can differ in the last digits between platforms.
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(script, out) {
  if (system2(rscript, script, stdout = out) != 0) stop(script, " failed", call. = FALSE)
}
same_csv <- function(f, tol) isTRUE(all.equal(read.csv(file.path("out/results", f)),
                                              read.csv(file.path("results", f)), tolerance = tol))
unlink("out", recursive = TRUE); dir.create("out/results", recursive = TRUE)
setwd("out"); invisible(file.copy("../code", ".", recursive = TRUE))
run("code/uls_identity.R", "results/uls_identity.txt")
if ("--dwls" %in% commandArgs(TRUE)) run("code/dwls_weights.R", "results/dwls_weights.txt")
setwd("..")

m <- read.csv("out/results/uls_models.csv")
ok <- c("uls_models.csv"                 = same_csv("uls_models.csv", 1e-6),
        "table1_rmsea_equivalent.csv"    = same_csv("table1_rmsea_equivalent.csv", 1e-10),
        "SRMR = sqrt(G/D)"               = max(abs(m$srmr_lavaan - m$srmr_formula)) < 1e-8,
        "SRMR = RMSEA sqrt(df/D)"        = max(abs(m$srmr_formula - m$rmsea_uls_formula * sqrt(m$df / m$D))) < 1e-12,
        "no model with SRMR > RMSEA"     = all(m$srmr_formula <= m$rmsea_uls_formula),
        "CFI = 1 - G/B (N = 1e5)"        = max(abs(m$cfi_lavaan - m$cfi_formula)) < 1e-3,
        "RMSEA, TLI never most lenient"  = all(m$most_lenient_clause %in% c("srmr", "cfi")))
if ("--dwls" %in% commandArgs(TRUE)) {
  d <- read.csv("out/results/dwls_weights.csv")
  ok <- c(ok, "dwls_weights.csv" = same_csv("dwls_weights.csv", 1e-4),
          "chisq = (N-1) F" = max(abs(d$chisq_over_N1F - 1)) < 1e-8,
          "RMSEA identity"  = max(abs(d$rmsea - d$rmsea_from_identity)) < 1e-8,
          "bounds hold"     = all(d$lower <= d$sqrt_F_over_df & d$sqrt_F_over_df <= d$upper))
}
for (f in names(ok)) cat(sprintf("%-32s %s\n", f, if (ok[[f]]) "reproduced" else "DIFFERS"))
if (!all(ok)) quit(status = 1)
