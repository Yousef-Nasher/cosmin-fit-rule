# population.R -- population values for the 38 conditions (no sampling).
# Run from study2/code/:  Rscript population.R   (writes out/population.csv)
# Fit indices under ULS and ML are computed from their definitions and checked against
# lavaan. Severity is measured on the parameters a user of the scale acts on (omega,
# factor correlations, loadings), never on fit indices, so the evaluation of the rule is
# not circular. A misspecification is "consequential" if omega is overstated by >= .05, a
# factor correlation or a loading is biased by >= .10, or a dimension correlated <= .70
# with its parent factor is omitted.
suppressPackageStartupMessages(library(lavaan)); source("design.R")
C <- conditions(); out <- list()
for (c in C) {
  B <- build_pop(c$p, c$m, c$load, c$phi, c$kind, c$k, c$delta)
  S <- B$Sig; p <- c$p; D <- p*(p+1)/2; one <- rep(1, p)
  omega_true <- as.numeric(t(one) %*% B$common %*% one) / sum(S)
  G0 <- sum(S[lower.tri(S)]^2); df0 <- p*(p-1)/2
  for (est in c("ULS","ML")) {
    f <- cfa(model_syntax(p, c$m), sample.cov = S, sample.nobs = 1e5,
             estimator = est, std.lv = TRUE)
    Sh <- fitted(f)$cov[rownames(S), colnames(S)]
    R  <- S - Sh; G <- sum(R[lower.tri(R)]^2); df <- as.numeric(fitMeasures(f, "df"))
    if (est == "ULS") { FM <- G; F0 <- G0 } else {
      FM <- log(det(Sh)) - log(det(S)) + sum(diag(S %*% solve(Sh))) - p
      F0 <- -log(det(S)) }
    sd <- lavInspect(f, "std"); lam <- sd$lambda; psi <- sd$psi
    lam_home <- rowSums(lam * simple_pattern(p, c$m))
    omega_fit <- sum(lam %*% psi %*% t(lam)) / sum(Sh)
    out[[length(out)+1]] <- data.frame(cond = c$id, p = p, m = c$m, load = c$load,
      kind = c$kind, k = c$k, delta = c$delta, est = est, df = df,
      rmsea_pop = sqrt(max(FM, 0) / df), srmr_pop = sqrt(G / D),
      cfi_pop = 1 - FM / F0, tli_pop = 1 - (FM/df) / (F0/df0),
      rmsea_lavaan = as.numeric(fitMeasures(f, "rmsea")),
      srmr_lavaan  = as.numeric(fitMeasures(f, "srmr")),
      omega_true = omega_true, omega_fit = omega_fit,
      omega_over = omega_fit - omega_true,
      fcor_bias  = if (c$m > 1) psi[1,2] - B$Phi[1,2] else NA_real_,
      load_bias  = max(abs(lam_home - c$load)))
  }
}
P <- do.call(rbind, out)
P$consequential <- (P$omega_over >= .05) | (abs(P$fcor_bias) >= .10 & !is.na(P$fcor_bias)) |
                   (P$load_bias >= .10) | (P$kind == "dim" & P$delta <= .70)
dir.create("out", showWarnings = FALSE); write.csv(P, "out/population.csv", row.names = FALSE)
cat("identity checks (ULS): max |rmsea_pop - sqrt(G/df)| is exact by construction;\n")
U <- P[P$est == "ULS", ]
cat(sprintf("  max |srmr_pop - lavaan srmr| = %.2e   (confirms D = p(p+1)/2)\n",
            max(abs(U$srmr_pop - U$srmr_lavaan))))
cat(sprintf("  SRMR <= RMSEA at population (ULS): %d / %d conditions\n",
            sum(U$srmr_pop <= U$rmsea_pop + 1e-12), nrow(U)))
cat(sprintf("  RMSEA<.06 & SRMR>=.08 at population, ULS: %d | ML: %d\n",
            sum(U$rmsea_pop < .06 & U$srmr_pop >= .08),
            sum(P$est=="ML" & P$rmsea_pop < .06 & P$srmr_pop >= .08)))
cat(sprintf("  consequential conditions: ULS %d/%d, ML %d/%d\n\n",
            sum(U$consequential), nrow(U), sum(P$consequential[P$est=="ML"]), sum(P$est=="ML")))
options(width = 160)
print(U[, c("cond","df","rmsea_pop","srmr_pop","cfi_pop","omega_over","fcor_bias","load_bias","consequential")],
      digits = 3, row.names = FALSE)
