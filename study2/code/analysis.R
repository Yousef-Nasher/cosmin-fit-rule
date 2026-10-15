# analysis.R -- Study 2: main analysis.
# Run from study2/:  Rscript code/analysis.R > results/study2_results.txt
# Reads the unzipped simulation output (out/final/*.csv) and out/population.csv.
options(width = 150)
source("code/rules.R")
D <- read_study2()
reps <- tapply(D$rep, D$cell, function(x) length(unique(x)))
cat(sprintf("chunks %d | rows %d | cells %d | reps per cell: min %d max %d\n\n",
    attr(D, "n_files"), nrow(D), length(reps), min(reps), max(reps)))
# ---- fit quality -------------------------------------------------------
cat("=== Fit quality by arm and N (percent of attempted fits) ===\n")
q <- aggregate(cbind(ok = D$status == "ok", admiss = D$status == "ok" & D$admiss %in% TRUE,
                     degenerate = D$status == "degenerate", nonconv = D$status == "nonconv") ~ arm + n,
               data = D, FUN = function(x) round(100 * mean(x), 2))
print(q, row.names = FALSE)
E <- clean_fits(D)

summ <- function(E, version, label, arms = c("ml_cont", "ml_ord_as_cont", "wlsmv_ord")) {
  R <- rules(idx(E, version), E$df, E$Dd)
  cat(sprintf("\n==================== %s ====================\n", label))
  for (a in arms) {
    cat(sprintf("\n--- arm: %s ---\n", a))
    cat(sprintf("%5s | %-26s | %-26s | %-26s | %s\n", "N", "COSMIN  false-fail/false-pass",
        "2-of-3  ff/fp", "all-3   ff/fp", "SRMR-df-scaled ff/fp"))
    for (n in sort(unique(E$n))) {
      g0 <- E$arm == a & E$n == n & E$kind == "none"
      g1 <- E$arm == a & E$n == n & E$conseq
      cell <- function(r) sprintf("%5.1f(%.1f) / %5.1f(%.1f)",
        100*mean(!R[[r]][g0]), 100*mcse(!R[[r]][g0]), 100*mean(R[[r]][g1]), 100*mcse(R[[r]][g1]))
      cat(sprintf("%5d | %-26s | %-26s | %-26s | %s\n", n, cell("COSMIN"),
          cell("two_of_three"), cell("all_three"), cell("SRMR_df_scaled")))
    }
    cat("  mechanism among ALL clean fits (percent):\n")
    for (n in sort(unique(E$n))) {
      g <- E$arm == a & E$n == n
      gp <- g & R$COSMIN & E$conseq     # false passes
      cat(sprintf("   N=%4d  RMSEA-rescues %5.2f | SRMR>RMSEA %5.1f | among false passes: SRMR sole %5.1f, CFI/TLI sole %5.1f | dropping RMSEA changes %d of %d ratings\n",
        n, 100*mean(R$rmsea_rescue[g]), 100*mean(R$srmr_gt_rmsea[g]),
        100*mean(R$srmr_sole[gp]), 100*mean(R$cfitli_sole[gp]),
        sum(R$COSMIN[g] != R$COSMIN_noRMSEA[g]), sum(g)))
    }
  }
  invisible(R)
}
R_sc  <- summ(E, "scaled", "PRIMARY: WLSMV scaled indices")
R_rob <- summ(E, "robust", "SENSITIVITY: WLSMV robust indices", arms = "wlsmv_ord")

cat("\n=== Which RMSEA? WLSMV scaled vs robust disagree on RMSEA<.06 (percent of fits) ===\n")
W <- E[E$arm == "wlsmv_ord" & !is.na(E$rmsea.rob), ]
dis <- (W$rmsea.sc < .06) != (W$rmsea.rob < .06)
for (n in sort(unique(W$n))) {
  cat(sprintf("  N=%4d  correct models: %5.1f   consequential: %5.1f\n", n,
      100*mean(dis[W$n == n & W$kind == "none"]), 100*mean(dis[W$n == n & W$conseq]))) }
cat(sprintf("  (robust RMSEA missing for %d of %d WLSMV fits)\n",
    sum(E$arm == "wlsmv_ord" & is.na(E$rmsea.rob)), sum(E$arm == "wlsmv_ord")))

cat("\n=== Rule results split by scale length (COSMIN false-pass / SRMR-df-scaled false-pass), N pooled ===\n")
for (a in c("ml_cont", "ml_ord_as_cont", "wlsmv_ord")) for (pp in sort(unique(E$p))) {
  g <- E$arm == a & E$p == pp & E$conseq
  if (sum(g) == 0) next
  cat(sprintf("  %-15s p=%2d  n=%5d  COSMIN %5.1f  SRMRdf %5.1f  2of3 %5.1f\n", a, pp, sum(g),
      100*mean(R_sc$COSMIN[g]), 100*mean(R_sc$SRMR_df_scaled[g]), 100*mean(R_sc$two_of_three[g]))) }

cat("\n=== Dose-response: COSMIN pass rate by population omega overstatement (resid conditions, all N) ===\n")
Rr <- R_sc; br <- c(-Inf, .02, .04, .06, .08, .12, Inf)
lab <- c("<.02",".02-.04",".04-.06",".06-.08",".08-.12",">.12")
E$ob <- cut(E$omega_over, br, labels = lab)
for (a in c("ml_cont", "ml_ord_as_cont", "wlsmv_ord")) {
  g <- E$arm == a & E$kind == "resid"
  cat(sprintf("  %-15s", a))
  for (b in lab) { h <- g & E$ob == b
    cat(sprintf(" %s:%5.1f/%5.1f", b, 100*mean(Rr$COSMIN[h]), 100*mean(Rr$all_three[h]))) }
  cat("   (COSMIN/all-3)\n") }

cat("\n=== Severity sensitivity: COSMIN false-pass with ML-based severity flag ===\n")
for (a in c("ml_cont", "ml_ord_as_cont", "wlsmv_ord")) for (n in sort(unique(E$n))) {
  g <- E$arm == a & E$n == n & E$conseq_ml
  cat(sprintf("  %-15s N=%4d  %5.1f\n", a, n, 100*mean(R_sc$COSMIN[g]))) }
saveRDS(list(E = E, R_sc = R_sc, R_rob = R_rob), "out/study2_analysis.rds")
