# dwls_identity.R -- quick test of the relations of study1/code/dwls_weights.R on small samples.
# Run from the repository root:  Rscript tests/dwls_identity.R      (needs lavaan; a few seconds)
# In any sample fitted by DWLS with free thresholds:
#   chisq = (N - 1) sum_i w_i e_i^2,   t SRMR^2 = sum of squared correlation residuals,
#   RMSEA^2 = w_eff (t / df) SRMR^2 - 1 / (N - 1),
#   sqrt(w_min t / df) SRMR <= sqrt(F / df) <= sqrt(w_max t / df) SRMR.
suppressPackageStartupMessages(library(lavaan))
check <- function(cond, label) {
  if (!isTRUE(cond)) stop("FAILED: ", label, call. = FALSE)
  cat("ok   ", label, "\n")
}
set.seed(20261006)
for (cuts in list(c(0), c(0.4, 1.5), c(-1.5, -0.5, 0.5, 1.5))) for (N in c(300, 1000)) {
  p <- 6; nm <- paste0("x", 1:p)
  S <- matrix(0.45^2, p, p); S[1, 2] <- S[2, 1] <- S[1, 2] + 0.30; diag(S) <- 1
  Z <- matrix(rnorm(N * p), N, p) %*% chol(S)
  X <- as.data.frame(lapply(seq_len(p), function(j) ordered(findInterval(Z[, j], cuts)))); names(X) <- nm
  fit <- cfa(paste("F1 =~", paste(nm, collapse = " + ")), data = X, ordered = nm,
             estimator = "WLSMV", std.lv = TRUE)
  e <- lavInspect(fit, "wls.obs") - lavInspect(fit, "wls.est")
  w <- diag(lavInspect(fit, "wls.v"))
  icor <- (length(e) - p * (p - 1) / 2 + 1):length(e)
  fm <- fitMeasures(fit, c("chisq", "df", "rmsea", "srmr"))
  t <- p * (p + 1) / 2; df <- fm[["df"]]; F <- sum(w * e^2); G <- sum(e[icor]^2)
  weff <- sum(w[icor] * e[icor]^2) / G; srmr <- sqrt(G / t)
  tag <- sprintf("(%d categories, N = %d)", length(cuts) + 1, N)
  check(max(abs(e[-icor])) < 1e-6, paste("threshold residuals are zero", tag))
  check(abs(fm[["chisq"]] - (N - 1) * F) < 1e-6 * fm[["chisq"]], paste("chisq = (N - 1) F", tag))
  check(abs(fm[["srmr"]] - srmr) < 1e-8, paste("SRMR from the correlation residuals", tag))
  check(abs(fm[["rmsea"]] - sqrt(max(weff * t * srmr^2 / df - 1 / (N - 1), 0))) < 1e-6,
        paste("RMSEA from w_eff and SRMR", tag))
  check(sqrt(min(w[icor]) * t / df) * srmr <= sqrt(F / df) + 1e-12 &&
        sqrt(F / df) <= sqrt(max(w[icor]) * t / df) * srmr + 1e-12, paste("bounds hold", tag))
}
cat("all DWLS identity checks passed\n")
