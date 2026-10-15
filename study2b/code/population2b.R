# population2b.R -- Study 2b: population fit and severity on the latent-response scale.
# Run from study2b/code/:  Rscript population2b.R   (writes out/population2b.csv; the frozen zip holds a copy)
# These values do not depend on the thresholds, so they are computed once per scale length,
# loading and model condition. Severity rule as in Study 2.
suppressPackageStartupMessages(library(lavaan)); source("design2b.R")
out <- list()
for (c in conditions()) {
  if (!(c$ct %in% c("c5_sym"))) next           # latent values do not depend on thresholds
  B <- build_pop(c$p, c$load, c$kind, c$k, c$delta); S <- B$Sig; p <- c$p; one <- rep(1, p)
  omega_true <- as.numeric(t(one) %*% B$common %*% one) / sum(S)
  for (est in c("ULS", "ML")) {
    f <- cfa(model_syntax(p), sample.cov = S, sample.nobs = 1e5, estimator = est, std.lv = TRUE)
    Sh <- fitted(f)$cov[rownames(S), colnames(S)]; R <- S - Sh
    G <- sum(R[lower.tri(R)]^2); df <- as.numeric(fitMeasures(f, "df"))
    FM <- if (est == "ULS") G else log(det(Sh)) - log(det(S)) + sum(diag(S %*% solve(Sh))) - p
    sd <- lavInspect(f, "std"); lam <- sd$lambda
    omega_fit <- sum(lam %*% t(lam)) / sum(Sh)
    out[[length(out) + 1]] <- data.frame(p = p, load = c$load, kind = c$kind, est = est, df = df,
      rmsea_pop = sqrt(max(FM, 0) / df), srmr_pop = sqrt(G / (p * (p + 1) / 2)),
      omega_over = omega_fit - omega_true, load_bias = max(abs(lam - c$load)))
  }
}
P <- do.call(rbind, out)
P$consequential <- (P$omega_over >= .05) | (P$load_bias >= .10) | (P$kind == "dim")
dir.create("out", showWarnings = FALSE); write.csv(P, "out/population2b.csv", row.names = FALSE)
print(P[P$est == "ULS", ], digits = 3, row.names = FALSE)
