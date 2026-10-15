# uls_identity.R -- Study 1, part A: population fit indices under unweighted least squares.
# Run from study1/:  Rscript code/uls_identity.R > results/uls_identity.txt
#
# A CFA is fitted by ULS to a population correlation matrix. With
#   G    = sum of squared residual correlations
#   B    = sum of squared correlations below the diagonal (misfit of the baseline model)
#   D    = p(p+1)/2, df = model degrees of freedom, df_B = p(p-1)/2
# the population indices are
#   SRMR  = sqrt(G / D)                 (lavaan; Shi, Maydeu-Olivares & Rosseel 2020)
#   RMSEA = sqrt(G / df)
#   CFI   = 1 - G / B
#   TLI   = 1 - (G / df) / (B / df_B)
# so the four clauses of the COSMIN rule are four thresholds on the same quantity G:
#   SRMR < .08  <=>  G < .0064 D          RMSEA < .06  <=>  G < .0036 df
#   CFI  > .95  <=>  G < .05 B            TLI   > .95  <=>  G < .05 B df / df_B
# Because df < D and df < df_B, the RMSEA and TLI clauses can never be the most lenient.
# Needs lavaan and the local package cosminfit (R CMD INSTALL ../cosminfit).
suppressPackageStartupMessages({ library(lavaan); library(cosminfit) })
options(width = 150)
cat(sprintf("R %s | lavaan %s | cosminfit %s\n\n", getRversion(),
            packageVersion("lavaan"), packageVersion("cosminfit")))

# ---- (a) grid of misspecified population models --------------------------
# Each row: true structure (m_true factors, loadings, factor correlation,
# residual correlations, cross-loadings) and the fitted simple-structure CFA (m_fit).
g <- function(...) expand.grid(..., stringsAsFactors = FALSE)
grid <- rbind(
  # A: true two-factor model fitted as one factor
  cbind(family = "A_2f_as_1f", g(p = c(4, 6, 8, 12, 20, 30), fcor = c(.5, .7, .9),
        load = c(.5, .7)), m_true = 2, m_fit = 1, k_rc = 0, rc = 0, n_xl = 0),
  # B: one-factor model with k unmodelled residual correlations
  cbind(family = "B_1f_resid", rbind(
        g(p = c(4, 5, 6, 8, 10, 12, 15, 20, 30), k_rc = 1, rc = .3, load = .7),
        g(p = c(4, 5, 6, 8, 10, 12, 15, 20, 30), k_rc = 3, rc = .4, load = .7),
        g(p = c(4, 5, 6, 8, 10, 12, 15, 20, 30), k_rc = 2, rc = .2, load = .5)),
        m_true = 1, m_fit = 1, fcor = 0, n_xl = 0),
  # C: m-factor model with 1-2 unmodelled cross-loadings of .3 (factor corr .3)
  cbind(family = "C_mf_xload", rbind(
        g(m_true = 2, p = c(6, 8, 10, 12, 15, 20, 30), n_xl = 1:2),
        g(m_true = 3, p = c(12, 15, 20, 30), n_xl = 1:2),
        g(m_true = 4, p = c(12, 20, 30), n_xl = 1:2)),
        fcor = .3, load = .7, k_rc = 0, rc = 0, m_fit = NA),
  # D: m-factor model with 2 unmodelled residual correlations of .3 (factor corr .3)
  cbind(family = "D_mf_resid", rbind(
        g(m_true = 2, p = c(6, 10, 20, 30)), g(m_true = 3, p = c(12, 15, 30)),
        g(m_true = 4, p = c(12, 20, 30))),
        fcor = .3, load = .5, k_rc = 2, rc = .3, n_xl = 0, m_fit = NA)
)
grid$m_fit[is.na(grid$m_fit)] <- grid$m_true[is.na(grid$m_fit)]   # C, D: fit the true m
grid <- grid[, c("family", "p", "m_true", "m_fit", "load", "fcor", "k_rc", "rc", "n_xl")]

# Items split as evenly as possible over m factors (factor 1 gets the first items)
assign_items <- function(p, m) sort(rep_len(seq_len(m), p))

# Population correlation matrix (standardised items, unit variances)
make_sigma <- function(s) {
  p <- s$p; m <- s$m_true; fac <- assign_items(p, m)
  L <- matrix(0, p, m); L[cbind(1:p, fac)] <- s$load
  if (s$n_xl > 0) for (j in seq_len(s$n_xl)) {     # first item of factor j also loads on factor j+1
    item <- which(fac == j)[1]; L[item, (j %% m) + 1] <- .3
  }
  Phi <- matrix(s$fcor, m, m); diag(Phi) <- 1
  S <- L %*% Phi %*% t(L)
  theta <- 1 - diag(S)                               # unique variances give unit item variances
  Th <- diag(theta)
  pairs <- list(c(1, 2), c(p - 1, p), c(2, 3))       # pairs carrying residual correlations
  if (s$k_rc > 0) for (j in seq_len(s$k_rc)) {
    a <- pairs[[j]][1]; b <- pairs[[j]][2]
    Th[a, b] <- Th[b, a] <- s$rc * sqrt(theta[a] * theta[b])  # residual correlation = rc
  }
  S <- S + Th
  stopifnot(all(theta > 0), min(eigen(S, only.values = TRUE)$values) > 0)
  dimnames(S) <- list(paste0("x", 1:p), paste0("x", 1:p))
  S
}

# Simple-structure CFA syntax: m_fit factors, each item on one factor
cfa_syntax <- function(p, m) {
  fac <- assign_items(p, m)
  paste(sapply(1:m, function(f) paste0("f", f, " =~ ",
        paste0("x", which(fac == f), collapse = " + "))), collapse = "\n")
}

# ---- (b) fit each model with ULS to Sigma and compute the quantities ------
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  s <- grid[i, ]; S <- make_sigma(s)
  fit <- cfa(cfa_syntax(s$p, s$m_fit), sample.cov = S, sample.nobs = 100000,
             estimator = "ULS", std.lv = TRUE)
  r  <- resid(fit, type = "raw")$cov                 # = residual correlations (Sigma has unit diag)
  lt <- r[lower.tri(r, diag = TRUE)]
  G  <- sum(lt^2); D <- s$p * (s$p + 1) / 2
  B  <- sum(S[lower.tri(S)]^2); dfB <- s$p * (s$p - 1) / 2
  fm <- fitMeasures(fit, c("df", "srmr", "rmsea", "fmin", "chisq", "cfi", "tli",
                           "baseline.chisq", "baseline.df"))
  data.frame(s, converged = lavInspect(fit, "converged"),
             max_abs_diag_resid = max(abs(diag(r))),
             G = G, D = D, df = unname(fm["df"]), df_cosminfit = cfa_df(s$p, s$m_fit),
             srmr_formula = sqrt(G / D), rmsea_uls_formula = sqrt(G / fm[["df"]]),
             srmr_lavaan = unname(fm["srmr"]), rmsea_lavaan = unname(fm["rmsea"]),
             fmin = unname(fm["fmin"]), chisq = unname(fm["chisq"]),
             B = B, df_B = dfB, cfi_formula = 1 - G / B,
             tli_formula = 1 - (G / fm[["df"]]) / (B / dfB),
             cfi_lavaan = unname(fm["cfi"]), tli_lavaan = unname(fm["tli"]),
             baseline_chisq = unname(fm["baseline.chisq"]))
}))
res$ratio_rmsea_srmr <- res$rmsea_uls_formula / res$srmr_formula
write.csv(res, "results/uls_models.csv", row.names = FALSE)

cat("=== (a) Grid ===\n")
cat(sprintf("models: %d | all converged: %s | p range %d-%d | fitted m: %s\n",
    nrow(res), all(res$converged), min(res$p), max(res$p),
    paste(sort(unique(res$m_fit)), collapse = ",")))
print(table(family = res$family, m_fit = res$m_fit))
cat(sprintf("max |diagonal residual|: %.3g\n", max(res$max_abs_diag_resid)))
cat(sprintf("lavaan df == cosminfit::cfa_df for all models: %s\n", all(res$df == res$df_cosminfit)))

cat("\n=== (b) Checks ===\n")
cat(sprintf("max |lavaan SRMR - sqrt(G/D)|                 : %.3g\n",
    max(abs(res$srmr_lavaan - res$srmr_formula))))
cat(sprintf("max |SRMR_formula - RMSEA_formula*sqrt(df/D)| : %.3g\n",
    max(abs(res$srmr_formula - res$rmsea_uls_formula * sqrt(res$df / res$D)))))
cat(sprintf("max |srmr_to_rmsea(SRMR, p, df) - RMSEA|      : %.3g  (cosminfit)\n",
    max(abs(srmr_to_rmsea(res$srmr_formula, res$p, res$df) - res$rmsea_uls_formula))))
cat(sprintf("models with SRMR > RMSEA (formula)             : %d of %d\n",
    sum(res$srmr_formula > res$rmsea_uls_formula), nrow(res)))
cat(sprintf("RMSEA/SRMR ratio range                         : %.3f to %.3f (= sqrt(D/df))\n",
    min(res$ratio_rmsea_srmr), max(res$ratio_rmsea_srmr)))
cat(sprintf("max |ratio - sqrt(D/df)|                       : %.3g\n",
    max(abs(res$ratio_rmsea_srmr - sqrt(res$D / res$df)))))
cat(sprintf("population SRMR range %.4f-%.4f | RMSEA_ULS range %.4f-%.4f\n",
    min(res$srmr_formula), max(res$srmr_formula),
    min(res$rmsea_uls_formula), max(res$rmsea_uls_formula)))
cat("\nlavaan's ULS fmin vs G:\n")
cat(sprintf("  2*fmin / G   : min %.8f  max %.8f\n", min(2 * res$fmin / res$G), max(2 * res$fmin / res$G)))
cat(sprintf("  chisq / ((N-1)*G), N = 1e5 : min %.8f  max %.8f\n",
    min(res$chisq / (99999 * res$G)), max(res$chisq / (99999 * res$G))))
N <- 100000
cat(sprintf("  max relative |2*fmin - G*(N-1)/N| / G : %.3g\n",
    max(abs(2 * res$fmin - res$G * (N - 1) / N) / res$G)))
cat(sprintf("  max relative |chisq - 2*N*fmin| / chisq : %.3g\n",
    max(abs(res$chisq - 2 * N * res$fmin) / res$chisq)))
cat(sprintf("  max |lavaan RMSEA - sqrt(G/df)|        : %.3g\n",
    max(abs(res$rmsea_lavaan - res$rmsea_uls_formula))))
cat(sprintf("  max |lavaan RMSEA - sqrt(max(chisq-df,0)/(df*(N-1)))| : %.3g\n",
    max(abs(res$rmsea_lavaan - sqrt(pmax(res$chisq - res$df, 0) / (res$df * (N - 1)))))))

cat("\nPer-family summary (formula values):\n")
print(do.call(rbind, lapply(split(res, res$family), function(d) data.frame(
  family = d$family[1], n = nrow(d), srmr_min = round(min(d$srmr_formula), 4),
  srmr_max = round(max(d$srmr_formula), 4), rmsea_min = round(min(d$rmsea_uls_formula), 4),
  rmsea_max = round(max(d$rmsea_uls_formula), 4),
  ratio_min = round(min(d$ratio_rmsea_srmr), 3), ratio_max = round(max(d$ratio_rmsea_srmr), 3),
  srmr_lt_.08_and_rmsea_ge_.06 = sum(d$srmr_formula < .08 & d$rmsea_uls_formula >= .06)))),
  row.names = FALSE)

# ---- (c) CFI and TLI share the numerator G; the rule as one threshold on G --
cat("\n=== (c) CFI, TLI and the rule as a single threshold on G ===\n")
cat(sprintf("max |baseline chisq / ((N-1) B) - 1|                       : %.3g\n",
    max(abs(res$baseline_chisq / ((N - 1) * res$B) - 1))))
cfi_n <- 1 - pmax((N - 1) * res$G - res$df, 0) / pmax((N - 1) * res$B - res$df_B, (N - 1) * res$G - res$df, 0)
tli_n <- 1 - ((N - 1) * res$G / res$df - 1) / ((N - 1) * res$B / res$df_B - 1)
cat(sprintf("max |lavaan CFI - CFI from G and B at N = 1e5|             : %.3g\n", max(abs(res$cfi_lavaan - cfi_n))))
cat(sprintf("max |lavaan TLI - TLI from G and B at N = 1e5|             : %.3g\n", max(abs(res$tli_lavaan - tli_n))))
cat(sprintf("max |lavaan CFI - (1 - G/B)|                               : %.3g\n", max(abs(res$cfi_lavaan - res$cfi_formula))))
cat(sprintf("max |lavaan TLI - (1 - (G/df)/(B/df_B))|                   : %.3g\n", max(abs(res$tli_lavaan - res$tli_formula))))
thr <- data.frame(srmr = .08^2 * res$D, rmsea = .06^2 * res$df, cfi = .05 * res$B,
                  tli = .05 * res$B * res$df / res$df_B)
lenient <- names(thr)[max.col(as.matrix(thr), ties.method = "first")]
pass <- data.frame(srmr = res$srmr_formula < .08, rmsea = res$rmsea_uls_formula < .06,
                   cfi = res$cfi_formula > .95, tli = res$tli_formula > .95)
cat(sprintf("clause pass == (G < its threshold) in all models           : %s\n",
    all(pass == (res$G < thr))))
cat(sprintf("rule passes == (G < largest threshold) in all models       : %s\n",
    all((rowSums(pass) > 0) == (res$G < do.call(pmax, thr)))))
cat("Most lenient clause (largest threshold on G), number of models:\n")
print(table(factor(lenient, levels = names(thr))))
cat("Most lenient clause by loading size:\n")
print(table(loading = res$load, clause = factor(lenient, levels = names(thr))))
cat(sprintf("models where RMSEA < .06 while SRMR >= .08                 : %d\n",
    sum(pass$rmsea & !pass$srmr)))
cat(sprintf("models where TLI > .95 while CFI <= .95                    : %d\n",
    sum(pass$tli & !pass$cfi)))
cat(sprintf("models passing the rule: %d of %d (SRMR only %d, CFI only %d, both %d)\n",
    sum(rowSums(pass) > 0), nrow(res), sum(pass$srmr & !pass$cfi), sum(pass$cfi & !pass$srmr),
    sum(pass$srmr & pass$cfi)))
res$most_lenient_clause <- lenient
write.csv(res, "results/uls_models.csv", row.names = FALSE)

# ---- (d) Table 1: RMSEA equivalent of SRMR .08, and SRMR needed for RMSEA .06
cat("\n=== (d) Table 1 ===\n")
pm <- data.frame(p = c(4, 5, 6, 8, 10, 12, 15, 20, 30, 8, 10, 12, 12, 15),
                 m = c(rep(1, 9), 2, 2, 2, 3, 3))
tab <- rbind(rmsea_equivalent_table(p = pm$p[pm$m == 1], m = 1),
             rmsea_equivalent_table(p = pm$p[pm$m == 2], m = 2),
             rmsea_equivalent_table(p = pm$p[pm$m == 3], m = 3))
tab[["srmr_needed_for_rmsea_.06"]] <- rmsea_to_srmr(.06, tab$p, tab$df)   # = .06*sqrt(df/D)
# sensitivity: off-diagonal convention D = p(p-1)/2
off <- rbind(rmsea_equivalent_table(p = pm$p[pm$m == 1], m = 1, convention = "offdiag"),
             rmsea_equivalent_table(p = pm$p[pm$m == 2], m = 2, convention = "offdiag"),
             rmsea_equivalent_table(p = pm$p[pm$m == 3], m = 3, convention = "offdiag"))
tab$D_offdiag <- off$D
tab$rmsea_equivalent_offdiag <- off$rmsea_equivalent
tab[["srmr_needed_for_rmsea_.06_offdiag"]] <- rmsea_to_srmr(.06, tab$p, tab$df, "offdiag")
write.csv(tab, "results/table1_rmsea_equivalent.csv", row.names = FALSE)
print(format(tab, digits = 4), row.names = FALSE)

# ---- (e) Holzinger-Swineford: which D does lavaan's SRMR use? ------------
cat("\n=== (e) Holzinger-Swineford 1939, 3-factor CFA, ML ===\n")
hs <- cfa("visual =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed =~ x7 + x8 + x9", data = HolzingerSwineford1939)
p <- 9
srmr_from <- function(r, D) sqrt(sum(r[lower.tri(r, diag = TRUE)]^2) / D)
r_bentler <- resid(hs, type = "cor.bentler")$cov   # residuals scaled by observed SDs
r_bollen  <- resid(hs, type = "cor.bollen")$cov    # observed minus implied correlations
cat(sprintf("lavaan srmr (default)              : %.6f\n", fitMeasures(hs, "srmr")))
cat(sprintf("lavaan srmr_bentler / srmr_mplus   : %.6f / %.6f\n",
    fitMeasures(hs, "srmr_bentler"), fitMeasures(hs, "srmr_mplus")))
cat(sprintf("max |diag| of cor.bentler / cor.bollen residuals: %.3g / %.3g\n",
    max(abs(diag(r_bentler))), max(abs(diag(r_bollen)))))
cat(sprintf("recomputed, cor.bentler, D = p(p+1)/2 = %d : %.6f\n", p * (p + 1) / 2, srmr_from(r_bentler, p * (p + 1) / 2)))
cat(sprintf("recomputed, cor.bentler, D = p(p-1)/2 = %d : %.6f\n", p * (p - 1) / 2, srmr_from(r_bentler, p * (p - 1) / 2)))
cat(sprintf("recomputed, cor.bollen,  D = p(p+1)/2 = %d : %.6f\n", p * (p + 1) / 2, srmr_from(r_bollen, p * (p + 1) / 2)))
cat(sprintf("recomputed, cor.bollen,  D = p(p-1)/2 = %d : %.6f\n", p * (p - 1) / 2, srmr_from(r_bollen, p * (p - 1) / 2)))
