# dwls_weights.R -- Study 1, part B: SRMR and RMSEA under diagonally weighted least squares.
# Run from study1/:  Rscript code/dwls_weights.R > results/dwls_weights.txt
#
# DWLS minimises F = sum_i w_i e_i^2, where e_i are the residuals of the thresholds and
# polychoric correlations and w_i = 1 / (N x sampling variance of statistic i). With free
# thresholds only the correlation residuals are non-zero, so with t = p(p+1)/2
#   chisq = (N - 1) F                (unscaled DWLS statistic, as in lavaan)
#   t SRMR^2 = sum_i e_i^2           (lavaan's SRMR)
#   F = w_eff t SRMR^2,  w_eff = sum_i w_i e_i^2 / sum_i e_i^2   (residual-weighted mean weight)
# These hold in every sample. They give
#   RMSEA^2 = w_eff (t / df) SRMR^2 - 1/(N - 1)     (unscaled RMSEA, when positive)
#   sqrt(w_min t / df) SRMR  <=  sqrt(F / df)  <=  sqrt(w_max t / df) SRMR
# With equal weights of 1 this is the ULS identity of part A. The mean-and-variance scaled
# statistic is chisq / c + shift, with scaling factor c and a shift that does not grow with
# N, so for the scaled RMSEA the weights are in effect w_i / c:
#   RMSEA_scaled^2 = (w_eff / c) (t / df) SRMR^2 + (shift - df) / (df (N - 1))
#
# The script checks these relations and tabulates the weights for the item formats,
# loadings and scale lengths of Study 2b, on one large sample per condition (N = 200,000),
# so that sample values are close to population values.
# Needs lavaan.
suppressPackageStartupMessages(library(lavaan))
options(width = 200)
cat(sprintf("R %s | lavaan %s\n\n", getRversion(), packageVersion("lavaan")))

N     <- 200000
SEED0 <- 20261006
CUTS  <- list(c2_sym = c(0),                    c2_skew = c(1.2),          # thresholds of Study 2b
              c3_sym = c(-0.75, 0.75),          c3_skew = c(0.4, 1.5),
              c5_sym = c(-1.5, -0.5, 0.5, 1.5), c5_skew = c(-0.2, 0.6, 1.2, 1.9))
COND  <- data.frame(p = c(5, 10, 10, 20), kind = c("resid", "resid", "dim", "dim"))
LOADS <- c(0.45, 0.70)

# Population correlation matrix of the latent responses (as in Study 2b):
# "resid": one factor plus a residual correlation of .30 between items 1 and 2;
# "dim":   two factors correlated .70, the second holding the last half of the items.
build_pop <- function(p, load, kind) {
  Lam <- matrix(load, p, 1); Phi <- matrix(1, 1, 1)
  if (kind == "dim") {
    idx <- (p - p %/% 2 + 1):p
    Lam <- cbind(Lam, 0); Lam[idx, 1] <- 0; Lam[idx, 2] <- load
    Phi <- matrix(c(1, .70, .70, 1), 2, 2)
  }
  Sig <- Lam %*% Phi %*% t(Lam)
  if (kind == "resid") Sig[1, 2] <- Sig[2, 1] <- Sig[1, 2] + .30
  diag(Sig) <- 1
  Sig
}

one <- function(ct, load, p, kind, seed) {
  set.seed(seed)
  nm <- paste0("x", 1:p)
  Z  <- matrix(rnorm(N * p), N, p) %*% chol(build_pop(p, load, kind))
  X  <- as.data.frame(lapply(seq_len(p), function(j) ordered(findInterval(Z[, j], CUTS[[ct]]))))
  names(X) <- nm
  fit <- cfa(paste("F1 =~", paste(nm, collapse = " + ")), data = X, ordered = nm,
             estimator = "WLSMV", std.lv = TRUE)
  e <- lavInspect(fit, "wls.obs") - lavInspect(fit, "wls.est")   # thresholds, then correlations
  w <- diag(lavInspect(fit, "wls.v"))
  ncor <- p * (p - 1) / 2; icor <- (length(e) - ncor + 1):length(e)
  ec <- e[icor]; wc <- w[icor]
  fm <- fitMeasures(fit, c("chisq", "df", "rmsea", "srmr", "chisq.scaled", "chisq.scaling.factor",
                           "rmsea.scaled", "rmsea.robust"))
  t <- p * (p + 1) / 2; df <- fm[["df"]]; G <- sum(ec^2); F <- sum(w * e^2)
  cc <- fm[["chisq.scaling.factor"]]; shift <- lavInspect(fit, "test")$scaled.shifted$shift.parameter
  srmr <- sqrt(G / t); weff <- sum(wc * ec^2) / G
  data.frame(format = ct, load = load, p = p, kind = kind, df = df, seed = seed,
    converged = lavInspect(fit, "converged"),
    max_abs_threshold_resid = max(abs(e[-icor])),
    chisq_over_N1F = fm[["chisq"]] / ((N - 1) * F),
    srmr = fm[["srmr"]], srmr_from_resid = srmr,
    w_min = min(wc), w_eff = weff, w_max = max(wc),
    rmsea = fm[["rmsea"]],
    rmsea_from_identity = sqrt(max(weff * t * srmr^2 / df - 1 / (N - 1), 0)),
    sqrt_F_over_df = sqrt(F / df),
    lower = sqrt(min(wc) * t / df) * srmr, upper = sqrt(max(wc) * t / df) * srmr,
    scaling_c = cc, shift = shift,
    rmsea_scaled = fm[["rmsea.scaled"]],
    rmsea_scaled_from_identity = sqrt(max(weff / cc * t * srmr^2 / df + (shift - df) / (df * (N - 1)), 0)),
    rmsea_scaled_without_shift = sqrt(weff / cc * t / df) * srmr,
    rmsea_robust = fm[["rmsea.robust"]])
}

grid <- merge(merge(data.frame(format = names(CUTS)), data.frame(load = LOADS)), COND)
grid$format <- factor(grid$format, levels = names(CUTS))
grid <- grid[order(grid$format, grid$load, grid$p, grid$kind), ]
grid$seed <- SEED0 + seq_len(nrow(grid))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  one(as.character(grid$format[i]), grid$load[i], grid$p[i], as.character(grid$kind[i]), grid$seed[i])))
res$t <- res$p * (res$p + 1) / 2
res$ratio_unscaled <- res$rmsea / res$srmr
res$ratio_scaled   <- res$rmsea_scaled / res$srmr
res$rmsea_scaled_at_srmr_.08 <- .08 * res$ratio_scaled
res$uls_ratio <- sqrt(res$t / res$df)
write.csv(res, "results/dwls_weights.csv", row.names = FALSE)

cat("=== (a) Checks, all conditions ===\n")
cat(sprintf("conditions: %d | all converged: %s | N = %d\n", nrow(res), all(res$converged), N))
cat(sprintf("max |threshold residual|                                   : %.3g\n", max(res$max_abs_threshold_resid)))
cat(sprintf("max |chisq / ((N-1) F) - 1|                                : %.3g\n", max(abs(res$chisq_over_N1F - 1))))
cat(sprintf("max |lavaan SRMR - sqrt(sum e^2 / t)|                      : %.3g\n", max(abs(res$srmr - res$srmr_from_resid))))
cat(sprintf("max |lavaan RMSEA - sqrt(w_eff (t/df) SRMR^2 - 1/(N-1))|   : %.3g\n", max(abs(res$rmsea - res$rmsea_from_identity))))
cat(sprintf("bound lower <= sqrt(F/df) <= upper in all conditions       : %s\n",
    all(res$lower <= res$sqrt_F_over_df & res$sqrt_F_over_df <= res$upper)))
cat(sprintf("max |lavaan scaled RMSEA - value from the identity|        : %.3g\n",
    max(abs(res$rmsea_scaled - res$rmsea_scaled_from_identity))))
cat(sprintf("max |lavaan scaled RMSEA - sqrt((w_eff/c) t/df) SRMR|      : %.3g   (shift term dropped)\n",
    max(abs(res$rmsea_scaled - res$rmsea_scaled_without_shift))))

cat("\n=== (b) Weights of the polychoric correlations, by item format and loading ===\n")
cat("(w = 1 / (N x sampling variance); ULS corresponds to w = 1; c = scaling factor of the scaled statistic)\n")
agg <- do.call(rbind, lapply(split(res, list(res$format, res$load), drop = TRUE), function(d) data.frame(
  format = d$format[1], load = d$load[1], w_min = min(d$w_min), w_eff_min = min(d$w_eff),
  w_eff_max = max(d$w_eff), w_max = max(d$w_max), c_min = min(d$scaling_c), c_max = max(d$scaling_c),
  w_eff_over_c_min = min(d$w_eff / d$scaling_c), w_eff_over_c_max = max(d$w_eff / d$scaling_c))))
agg <- agg[order(factor(agg$format, levels = names(CUTS)), agg$load), ]
print(format(agg, digits = 3), row.names = FALSE)

cat("\n=== (c) RMSEA relative to SRMR, by condition ===\n")
cat("(uls_ratio = sqrt(t/df), the RMSEA/SRMR ratio under ULS; the last column is the scaled RMSEA\n")
cat(" that corresponds to SRMR = .08 at this ratio)\n")
show <- res[, c("format", "load", "p", "kind", "df", "srmr", "rmsea", "rmsea_scaled", "rmsea_robust",
                "w_eff", "scaling_c", "uls_ratio", "ratio_unscaled", "ratio_scaled", "rmsea_scaled_at_srmr_.08")]
print(format(show, digits = 3), row.names = FALSE)

cat("\n=== (d) Where the RMSEA clause is the more lenient one ===\n")
cat("RMSEA < .06 is easier to meet than SRMR < .08 when RMSEA/SRMR < .75.\n")
tb <- function(x) table(format = res$format, load = res$load, x)
cat(sprintf("unscaled RMSEA below SRMR: %d of %d conditions; scaled RMSEA below SRMR: %d of %d\n",
    sum(res$ratio_unscaled < 1), nrow(res), sum(res$ratio_scaled < 1), nrow(res)))
cat(sprintf("scaled RMSEA / SRMR < .75: %d of %d conditions\n", sum(res$ratio_scaled < .75), nrow(res)))
print(res[res$ratio_scaled < 1, c("format", "load", "p", "kind", "ratio_unscaled", "ratio_scaled")], digits = 3, row.names = FALSE)
cat(sprintf("\nscaled RMSEA corresponding to SRMR = .08: range %.3f to %.3f over all conditions\n",
    min(res$rmsea_scaled_at_srmr_.08), max(res$rmsea_scaled_at_srmr_.08)))
print(round(tapply(res$rmsea_scaled_at_srmr_.08, list(format = res$format, paste0("p", res$p, "_", res$kind, "_l", res$load)), mean), 3))
