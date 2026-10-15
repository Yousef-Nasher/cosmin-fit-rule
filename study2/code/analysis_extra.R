# analysis_extra.R -- Study 2: what the RMSEA clause does, results per condition, and the
# sensitivity analysis that keeps inadmissible fits.
# Run from study2/ after analysis.R:  Rscript code/analysis_extra.R > results/study2_extra.txt
X <- readRDS("out/study2_analysis.rds"); E <- X$E; R <- X$R_sc; options(width = 200)
E$cond <- sub("_n[0-9]{4}$", "", E$cell)
# 1. What the RMSEA clause does: ratings that change when it is dropped, split by model truth
cat("=== Effect of the RMSEA clause (ratings changed by dropping it), WLSMV scaled primary ===\n")
chg <- R$COSMIN != R$COSMIN_noRMSEA
for (a in unique(E$arm)) for (n in 100) {
  g <- E$arm == a & E$n == n
  cat(sprintf("  %-15s N=%d  changed %4d of %5d | correct model rescued %4d | consequential passed %4d | non-consequential misfit passed %4d\n",
      a, n, sum(chg[g]), sum(g), sum(chg[g] & E$kind[g] == "none"), sum(chg[g] & E$conseq[g]),
      sum(chg[g] & E$kind[g] != "none" & !E$conseq[g])))
}
cat("\n  N=100, WLSMV: rescue rate of CORRECT models by loading and scale (percent)\n")
g <- E$arm == "wlsmv_ord" & E$n == 100 & E$kind == "none"
print(round(100 * tapply(chg[g], list(paste0("p", E$p[g], "_m", E$m[g]), load = E$load[g]), mean), 1))
cat("\n  same, robust indices (sensitivity)\n")
chg2 <- X$R_rob$COSMIN != X$R_rob$COSMIN_noRMSEA
print(round(100 * tapply(chg2[g], list(paste0("p", E$p[g], "_m", E$m[g]), load = E$load[g]), mean), 1))

# 2. Per-condition COSMIN pass rates (design-free view), WLSMV and ML continuous
cat("\n=== COSMIN pass rate per condition (percent), rows = condition, cols = N ===\n")
P <- read.csv("out/population.csv"); PU <- P[P$est == "ULS", ]
for (a in c("wlsmv_ord", "ml_cont")) {
  g <- E$arm == a
  tb <- round(100 * tapply(R$COSMIN[g], list(E$cond[g], E$n[g]), mean), 1)
  m <- match(rownames(tb), PU$cond)
  out <- data.frame(tb, conseq = PU$consequential[m], rmsea_pop = round(PU$rmsea_pop[m], 3),
                    srmr_pop = round(PU$srmr_pop[m], 3), omega_over = round(PU$omega_over[m], 3), check.names = FALSE)
  cat(sprintf("\n--- %s ---\n", a)); print(out)
}
# 3. Sensitivity: include inadmissible converged fits
cat("\n=== Sensitivity: false-fail / false-pass including inadmissible converged fits (WLSMV scaled) ===\n")
fs <- list.files("out/final", "\\.csv$", full.names = TRUE)
D <- do.call(rbind, lapply(fs, read.csv)); D$cond <- sub("_n[0-9]{4}$", "", D$cell)
D <- D[D$status == "ok" & D$arm == "wlsmv_ord", ]; D$conseq <- PU$consequential[match(D$cond, PU$cond)]
f <- function(v) !is.na(v) & v
pass <- f(D$cfi.sc > .95) | f(D$tli.sc > .95) | f(D$rmsea.sc < .06) | f(D$srmr < .08)
for (n in sort(unique(D$n))) cat(sprintf("  N=%4d  false-fail %.1f  false-pass %.1f   (inadmissible fits: %d)\n", n,
  100 * mean(!pass[D$n == n & D$kind == "none"]), 100 * mean(pass[D$n == n & D$conseq]), sum(D$n == n & !(D$admiss %in% TRUE))))
